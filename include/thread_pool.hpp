/**
 * @file thread_pool.hpp
 * @brief A simple thread pool which automatically execute tasks added, suspend
 * and stop.
 * @author Haoming Bai <haomingbai@hotmail.com>
 * @date   2025-06-19
 *
 * Copyright © 2025 Haoming Bai
 * SPDX-License-Identifier: MIT
 *
 * @details
 * This thread pool implementation provides:
 * 1. Fixed-size worker threads that process tasks from a thread-safe queue
 * 2. Three operational states: RUN (normal processing), SUSPEND (paused), STOP
 * (shutdown)
 * 3. Full lifecycle control with pause/resume/stop functionality
 * 4. Exception handling callback for task exceptions
 * 5. Future-based task submission with result retrieval
 *
 * Key synchronization features:
 * - Condition variable manages thread suspension/resumption
 * - Atomic state transitions ensure consistent pool status
 * - Worker threads block efficiently during suspension/empty queues
 * - Queue notifications minimize task processing latency
 *
 * The pool guarantees:
 * - Tasks posted via post() return futures for result retrieval
 * - Stopped pools join all threads upon destruction
 * - Exception handlers capture std::exception messages
 * - Custom queue types may be used if satisfying SafeQueueLike concept
 */

#pragma once
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>

#include "safe_queue.hpp"

namespace blib {

/**
 * @brief A simple thread pool implementation with pause, resume, and stop
 * controls.
 *
 * The ThreadPool class manages a fixed number of worker threads that
 * continuously fetch and execute tasks from an internal thread-safe queue.
 * Tasks are posted as function objects and executed asynchronously.
 */
template <typename Queue = SafeQueue<std::function<void()>>>
  requires SafeQueueLike<Queue, std::function<void()>>
class ThreadPool {
  Queue tsk_queue_;
  std::vector<std::thread> threads_;
  std::mutex mtx_;
  std::condition_variable cv_;
  enum status { STOP, SUSPEND, RUN };
  std::atomic<status> stat_;
  std::function<void(const char *)> handle_exception_;

 public:
  /**
   * @brief Construct a ThreadPool with given number of threads.
   *
   * @param size Number of worker threads to launch.
   * Threads start immediately and wait for tasks.
   */
  explicit ThreadPool(size_t size)
      : threads_(size), handle_exception_([](const char *) {}) {
    stat_ = RUN;
    for (auto &it : threads_) {
      it = std::thread([this] {
        while (true) {
          if (this->stat_.load() == STOP) {
            break;
          }
          if (this->stat_.load() == SUSPEND) {
            std::unique_lock<std::mutex> lck(this->mtx_);
            cv_.wait(lck, [this] { return stat_.load() != SUSPEND; });
            continue;
          }
          std::function<void()> fn;
          auto succ = tsk_queue_.pop(fn);
          if (succ) {
            try {
              fn();
            } catch (const std::exception &e) {
              handle_exception_(e.what());
            } catch (...) {
            }
          } else {
            std::unique_lock<std::mutex> lock(mtx_);
            cv_.wait(lock);
          }
        }
      });
    }
  }

  /**
   * @brief Construct a ThreadPool with exception handling callback.
   *
   * @param size Number of worker threads to launch.
   * @param exception_handler Callback invoked with exception message when a
   * task throws.
   */
  ThreadPool(size_t size, std::function<void(const char *)> exception_handler)
      : threads_(size), handle_exception_(exception_handler) {
    stat_ = RUN;
    for (auto &it : threads_) {
      it = std::thread([this] {
        while (true) {
          if (this->stat_.load() == STOP) {
            break;
          }
          if (this->stat_.load() == SUSPEND) {
            std::unique_lock<std::mutex> lck(this->mtx_);
            cv_.wait(lck, [this] { return stat_.load() != SUSPEND; });
            continue;
          }
          std::function<void()> fn;
          auto succ = tsk_queue_.pop(fn);
          if (succ) {
            try {
              fn();
            } catch (const std::exception &e) {
              handle_exception_(e.what());
            } catch (...) {
            }
          } else {
            std::unique_lock<std::mutex> lock(mtx_);
            cv_.wait(lock);
          }
        }
      });
    }
  }

  /**
   * @brief Stop all threads and prevent further task execution.
   *
   * Notifies all worker threads and joins them upon destruction or explicit
   * join().
   */
  void stop() {
    if (stat_.load() != STOP) {
      stat_ = STOP;
      cv_.notify_all();
    }
  }

  /**
   * @brief Resume or restart the thread pool.
   *
   * - If the pool was stopped, joins old threads and creates new ones.
   * - If the pool was suspended, wakes all threads to RUN state.
   */
  void run() {
    if (stat_.load() == STOP) {
      for (auto &it : threads_) {
        if (it.joinable()) {
          it.join();
        }
      }
      stat_ = RUN;
      for (auto &it : threads_) {
        it = std::thread([this] {
          while (true) {
            if (this->stat_.load() == STOP) {
              break;
            }
            if (stat_.load() == SUSPEND) {
              std::unique_lock<std::mutex> lck(this->mtx_);
              cv_.wait(lck, [this] { return stat_.load() != SUSPEND; });
            }
            std::function<void()> fn;
            auto succ = tsk_queue_.pop(fn);
            if (succ) {
              try {
                fn();
              } catch (const std::exception &e) {
                handle_exception_(e.what());
              } catch (...) {
              }
            } else {
              std::unique_lock<std::mutex> lock(mtx_);
              cv_.wait(lock);
            }
          }
        });
      }
    } else if (this->stat_.load() == SUSPEND) {
      stat_ = RUN;
      this->cv_.notify_all();
    }
  }

  /**
   * @brief Temporarily suspend task execution.
   *
   * Worker threads will block until run() is called.
   */
  void suspend() {
    if (stat_.load() == RUN) {
      this->stat_ = SUSPEND;
    }
  }

  /**
   * @brief Stop and join all worker threads.
   *
   * Equivalent to stop() followed by joining each thread.
   */
  void join() {
    this->stop();
    for (auto &it : threads_) {
      if (it.joinable()) {
        it.join();
      }
    }
  }

  /**
   * @brief post a task to the thread pool.
   *
   * Accepts a callable and its arguments, schedules it for
   * asynchronous execution, and returns a future for the result.
   *
   * @tparam Fn Callable type.
   * @tparam Args Argument types.
   * @param func Callable object or function pointer.
   * @param args Arguments to pass to the callable.
   * @return std::future holding the return value of the callable.
   */
  template <typename Fn, typename... Args>
  auto post(Fn &&func, Args &&...args) -> std::future<decltype(func(args...))> {
    std::function<decltype(func(args...))()> fn_without_param =
        std::bind(std::forward<Fn>(func), std::forward<Args...>(args)...);
    auto task_ptr =
        std::make_shared<std::packaged_task<decltype(func(args...))()>>(
            fn_without_param);
    std::function<void()> wrapper = [task_ptr] { (*task_ptr)(); };

    this->tsk_queue_.push(wrapper);
    cv_.notify_one();

    return task_ptr->get_future();
  }

  /**
   * @brief Destructor: stops and joins all threads.
   */
  ~ThreadPool() { this->join(); }

  /// Disable copy construction and move construction.
  ThreadPool(ThreadPool &&) = delete;
  ThreadPool(const ThreadPool &) = delete;
  ThreadPool &operator=(ThreadPool &&) = delete;
  ThreadPool &operator=(const ThreadPool &) = delete;
};

}  // namespace blib

#endif  // THREAD_POOL_HPP
