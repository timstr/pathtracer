#include <Renderer.hpp>

#include <RandomNumberGenerator.hpp>

#include <cassert>
#include <iostream>
#include <random>

Renderer::Renderer() noexcept
    : m_renderBarrierMaybe(std::nullopt)
    , m_renderData(std::nullopt) {
}

Renderer::~Renderer() {
    this->stopThreadPool();
}

void Renderer::startThreadPool(size_t numThreads) {
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
    }
    std::cout << "Starting thread pool using " << numThreads << " threads"
        << std::endl;
    this->stopThreadPool();
    assert(!m_renderBarrierMaybe.has_value());
    assert(!m_renderData.has_value());
    assert(m_threadPool.size() == 0);
    // NOTE: the current thread must also arrive at the barrier to allow
    // phase transitions
    m_renderBarrierMaybe.emplace(numThreads + 1);
    m_timeToExit.store(false);
    for (std::size_t i = 0; i < numThreads; ++i) {
        m_threadPool.emplace_back([this]{ this->doWork(); });
    }
}

void Renderer::stopThreadPool() {
    if (m_threadPool.size() == 0) {
        assert(!m_renderBarrierMaybe.has_value());
        assert(!m_renderData.has_value());
        return;
    }
    assert(m_renderBarrierMaybe.has_value());
    assert(m_timeToExit.load() == false);
    m_timeToExit.store(true);
    m_renderBarrierMaybe->arrive_and_wait();
    for (auto& t : m_threadPool) {
        t.join();
    }
    m_threadPool.clear();
    m_renderBarrierMaybe.reset();
}

void Renderer::doWork() const noexcept {
    while (true) {
        // Phase 1: waiting for work
        assert(m_renderBarrierMaybe.has_value());
        m_renderBarrierMaybe->arrive_and_wait();
        if (m_timeToExit.load()) {
            return;
        }

        // Phase 2: doing the work
        assert(m_renderData.has_value());
        const auto& camera = *m_renderData->camera;
        const auto& scene = *m_renderData->scene;
        const auto& settings = *m_renderData->settings;
        auto& img = *m_renderData->image;
        assert(settings.width() == img.width());
        assert(settings.height() == img.height());

        const auto maxWidth = static_cast<float>(settings.width() - 1);
        const auto maxHeight = static_cast<float>(settings.height() - 1);
        const auto sppInv = 1.0f / static_cast<float>(settings.samplesPerPixel());
        const auto deltaX = 1.0f / maxWidth;
        const auto deltaY = 1.0f / maxHeight;
        const auto jitterX = std::uniform_real_distribution<float>{ -0.5f * deltaX, 0.5f * deltaX };
        const auto jitterY = std::uniform_real_distribution<float>{ -0.5f * deltaY, 0.5f * deltaY };
        while (true) {
            const auto taskIndex = m_nextTaskIndex.fetch_add(1);
            auto taskMaybe = Renderer::makeTask(taskIndex, settings.width(), settings.height());
            if (!taskMaybe.has_value()) {
                break;
            }
            const auto& task = *taskMaybe;
            const auto y = task.y;
            const auto py = static_cast<float>(y) / maxHeight;
            for (std::size_t x = task.xStart; x < task.xEnd; ++x) {
                const auto px = static_cast<float>(x) / maxWidth;
                auto acc = Color{};
                for (std::size_t i = 0; i < settings.samplesPerPixel(); ++i){
                    const auto sx = px + jitterX(randomEngine());
                    const auto sy = py + jitterY(randomEngine());
                    const auto viewRay = camera.getViewRay(sx, sy);
                    const auto c = scene.trace(viewRay, settings.numBounces());
                    // TODO: operators
                    acc.r += c.r;
                    acc.g += c.g;
                    acc.b += c.b;
                }
                // TODO: operators
                auto& pixel = img(x, y);
                pixel.r = acc.r * sppInv;
                pixel.g = acc.g * sppInv;
                pixel.b = acc.b * sppInv;
            }
        }

        assert(m_renderBarrierMaybe.has_value());
        m_renderBarrierMaybe->arrive_and_wait();
        if (m_timeToExit.load()) {
            return;
        }
    }
}

std::optional<Renderer::render_task> Renderer::makeTask(
    size_t taskIndex,
    size_t width,
    size_t height
) noexcept {
    const auto tasksPerRow = width / s_pixelsPerTask;
    const auto xStart = s_pixelsPerTask * (taskIndex % tasksPerRow);
    const auto y = taskIndex / tasksPerRow;
    if (y >= height) {
        return std::nullopt;
    }
    const auto xEnd = std::min(
        xStart + s_pixelsPerTask,
        width
    );
    return render_task {
        xStart,
        xEnd,
        y
    };
}

Image Renderer::render(
    const Scene& scene,
    const Camera& camera,
    const RenderSettings& settings
) {
    const auto lock = std::lock_guard(m_mutex);
    if (m_threadPool.size() == 0) {
        this->startThreadPool();
    }

    auto img = Image(settings.width(), settings.height());

    assert(!m_renderData.has_value());
    m_renderData.emplace(RenderData{
        &scene,
        &camera,
        &settings,
        &img
    });
    m_nextTaskIndex.store(0);
    assert(m_renderBarrierMaybe.has_value());
    assert(!m_timeToExit.load());

    // Phase 1: threads are all waiting at the barrier to start working.
    // This last thread arriving starts them
    m_renderBarrierMaybe->arrive_and_wait();
    // Phase 2: threads are all busy doing work.
    // When the last thread arrives, all work is done.
    m_renderBarrierMaybe->arrive_and_wait();

    m_renderData.reset();

    return img;
}

void Renderer::cancelRender() noexcept {
    m_nextTaskIndex.store(std::numeric_limits<size_t>::max() - m_threadPool.size());
}