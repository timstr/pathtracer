#pragma once

#include <Camera.hpp>
#include <Image.hpp>
#include <RenderSettings.hpp>
#include <Scene.hpp>

#include <atomic>
#include <barrier>
#include <cstddef>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

class Renderer {
public:
    Renderer() noexcept;
    ~Renderer();

    Renderer(const Renderer&) = delete;
    Renderer(Renderer&&) = delete;
    Renderer& operator=(const Renderer&) = delete;
    Renderer& operator=(Renderer&&) = delete;

    void startThreadPool(size_t numThreads = 0);
    void stopThreadPool();

    Image render(const Scene&, const Camera&, const RenderSettings&);

private:
    mutable std::mutex m_mutex;
    mutable std::optional<std::barrier<>> m_renderBarrierMaybe;
    mutable std::atomic<size_t> m_nextTaskIndex;
    std::atomic<bool> m_timeToExit;

    static const size_t s_pixelsPerTask = 8;

    std::vector<std::thread> m_threadPool;

    struct render_task {
        size_t xStart;
        size_t xEnd;
        size_t y;
    };

    static std::optional<render_task> makeTask(
        size_t taskIndex,
        size_t width,
        size_t height
    ) noexcept;

    struct RenderData {
        const Scene* scene;
        const Camera* camera;
        const RenderSettings* settings;
        Image* image;
    };

    std::optional<RenderData> m_renderData;

    void doWork() const noexcept;
};
