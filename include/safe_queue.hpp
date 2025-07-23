/**
 * @file safe_queue.hpp
 * @brief A thread-safe queue which can support multiple methods to fetch.
 * elements.
 * @author Haoming Bai <haomingbai@hotmail.com>
 * @date   2025-06-19
 *
 * Copyright © 2025 Haoming Bai
 * SPDX-License-Identifier: MIT
 *
 * @details
 * This thread-safe queue implementation guarantees safe concurrent access
 * using:
 * 1. A counting semaphore to block pop operations when empty
 * 2. A mutex-protected std::queue for underlying storage
 *
 * Key synchronization behavior:
 * - Blocking pop() suspends threads without CPU consumption
 * - push() always releases the semaphore after mutex unlock
 * - Non-blocking/timed pop variants avoid indefinite waits
 *
 * The semaphore tracks queued items independently of the actual queue size.
 * All operations maintain strong exception safety where possible.
 */

#pragma once
#ifndef SAFE_QUEUE_HPP
#define SAFE_QUEUE_HPP

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstddef>
#include <limits>
#include <mutex>
#include <queue>
#include <semaphore>
#include <utility>

namespace blib {

template <typename Q, typename T>
concept SafeQueueLike = requires(Q &q, T &&in, T &out) {
  { q.push(std::forward<T>(in)) };
  { q.pop(out) } -> std::convertible_to<bool>;
  { q.empty() } -> std::convertible_to<bool>;
};

/**
 * @brief A thread-safe queue that can block poprs when empty.
 *
 * SafeQueue uses a std::counting_semaphore to make `pop()` block
 * until an item is available, and a mutex to protect the underlying
 * std::queue.
 *
 * Unblocked methods are also available if put the reference as a param in the
 * pop() function.
 *
 * @tparam T Type of the elements stored in the queue.
 */
template <typename T>
class SafeQueue {
  std::queue<T> queue_;  ///< Underlying FIFO container, the queue should be
                         ///< decleared first to be safely constructed.
  std::counting_semaphore<std::numeric_limits<int>::max()>
      sem_;         ///< Counts available items
  std::mutex mtx_;  ///< Protects access to queue

 public:
  /**
   * @brief Construct an empty SafeQueue.
   *
   * The semaphore is initialized to zero, so `pop()` will block
   * until an element is pushd.
   */
  SafeQueue() : queue_(), sem_(0) {}

  /**
   * @brief Construct a SafeQueue from an existing queue.
   * @param oldQueue A pre-filled standard queue; contents are moved.
   *
   * The semaphore is initialized to the size of the moved-in queue,
   * so poprs will be unblocked that many times.
   */
  SafeQueue(std::queue<T> &&old_queue)
      : queue_(std::move(old_queue)), sem_(queue_.size()) {}

  // Disable copying and moving
  SafeQueue(SafeQueue &&) = delete;
  SafeQueue(const SafeQueue &) = delete;
  SafeQueue &operator=(SafeQueue &&) = delete;
  SafeQueue &operator=(const SafeQueue &) = delete;
  ~SafeQueue() = default;

  /**
   * @brief Remove and return the front element.
   *
   * If the queue is empty, blocks the calling thread until
   * another thread calls push().
   *
   * Notice: once the thread is blocked, it will not occupy the cpu but suspend
   * and wait aside until it is woked.
   *
   * @return The front element.
   */
  T pop() {
    sem_.acquire();
    std::lock_guard<std::mutex> lock(mtx_);
    auto ret = std::move(queue_.front());
    queue_.pop();
    return ret;
  }

  /**
   * @brief Try to pop an element without blocking.
   *
   * Attempts to acquire the semaphore immediately. If an item is available,
   * locks the queue, moves the front element into `obj`, pops it, and returns
   * true. Otherwise returns false immediately.
   *
   * @param[out] obj  Reference into which the popd element will be moved.
   * @return true if an element was popd, false if the queue was empty.
   */
  bool pop(T &obj) {
    auto res = sem_.try_acquire();
    if (res) {
      std::lock_guard<std::mutex> lock(mtx_);
      obj = std::move(queue_.front());
      queue_.pop();
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Try to pop an element, blocking up to a timeout.
   *
   * Attempts to acquire the semaphore, waiting up to the specified duration.
   * If an item becomes available within `time`, locks the queue, moves the
   * front element into `obj`, pops it, and returns true. If the timeout
   * expires without an item becoming available, returns false.
   *
   * @tparam Time   A std::ratio-based duration rep (e.g., std::milli).
   * @param[out] obj  Reference into which the popd element will be moved.
   * @param time      Maximum duration to wait for an available element.
   * @return true if an element was popd within the timeout, false if timed
   * out.
   */
  template <typename Time>
  bool pop(T &obj, std::chrono::duration<Time> time) {
    auto res = sem_.try_acquire_for(time);
    if (res) {
      std::lock_guard<std::mutex> lock(mtx_);
      obj = std::move(queue_.front());
      queue_.pop();
      return true;
    } else {
      return false;
    }
  }

  /**
   * @brief Add a new element to the back of the queue.
   *
   * This will un-block one thread waiting in pop().
   *
   * @tparam Args Types that can be forwarded to construct T.
   * @param args Arguments forwarded to T’s constructor.
   */
  template <typename... Args>
  void push(Args &&...args) {
    {
      std::lock_guard<std::mutex> lock(mtx_);
      queue_.emplace(std::forward<Args...>(args)...);
    }
    sem_.release();
  }

  /**
   * @brief Test whether the queue is empty.
   *
   * Non-blocking; just returns the current size == 0.
   *
   * @return true if empty, false otherwise.
   */
  bool empty() {
    std::lock_guard<std::mutex> lock(mtx_);
    return queue_.empty();
  }

  /**
   * @brief Get the number of elements pending.
   * @return Current queue size.
   */
  std::size_t size() {
    std::lock_guard<std::mutex> lock(mtx_);
    return queue_.size();
  }
};

}  // namespace blib

#endif  // SAFE_QUEUE_HPP
