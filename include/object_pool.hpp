/**
 * @file object_pool.hpp
 * @brief A generic thread-safe object pool template
 * @author Haoming Bai <haomingbai@hotmail.com>
 * @date   2025-06-25
 *
 * Copyright © 2025 Haoming Bai
 * SPDX-License-Identifier: MIT
 *
 * @details This header defines ObjectPool, a templated thread-safe pool that
 * manages the lifecycle of objects of type DataType. Clients can acquire and
 * release objects efficiently. The pool supports custom health checks, factory
 * functions, dynamic resizing, and optionally integrates with SafeQueue-like
 * queue implementations.
 */

#pragma once
#ifndef OBJECT_POOL_HPP
#define OBJECT_POOL_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>

#include "safe_queue.hpp"

namespace blib {

/**
 * @class ObjectPool
 * @brief Manages a pool of reusable objects in a thread-safe manner.
 * @tparam DataType Type of objects managed by the pool.
 * @tparam Queue Underlying queue type (must satisfy SafeQueueLike).
 *
 * The pool pre-allocates objects via a factory function and performs health
 * checks before issuing and upon return. Objects are wrapped in ObjectWrapper
 * to ensure automatic return to the pool when they go out of scope.
 */
template <typename DataType, typename Queue = SafeQueue<DataType>>
  requires SafeQueueLike<Queue, std::unique_ptr<DataType>>
class ObjectPool
    : public std::enable_shared_from_this<ObjectPool<DataType, Queue>> {
  std::function<bool(const DataType &)>
      check_health_;                  /**< Health check predicate. */
  std::function<DataType()> factory_; /**< Factory to create new instances. */
  Queue queue_;                       /**< Underlying thread-safe queue. */
  std::size_t size_; /**< Maximum pool size (initial capacity). */

 public:
  /**
   * @brief Constructs an empty pool with a health check and custom factory.
   * @param check_health Function to validate objects.
   * @param get_new Factory function to create new DataType instances.
   */
  ObjectPool(std::function<bool(const DataType &)> check_health,
             std::function<DataType()> get_new)
      : check_health_(check_health), factory_(get_new), size_(0) {}

  /**
   * @brief Constructs an empty pool with custom health check and default
   * factory.
   * @param check_health Function to validate objects.
   */
  ObjectPool(std::function<bool(const DataType &)> check_health)
      : check_health_(check_health), factory_([] { return DataType{}; }) {}

  /**
   * @brief Constructs an empty pool with default health check and custom
   * factory.
   * @param get_new Factory function to create new DataType instances.
   */
  explicit ObjectPool(std::function<DataType()> get_new)
      : check_health_([](const DataType &obj) { return true; }),
        factory_(get_new),
        size_(0) {}

  /**
   * @brief Constructs an empty pool with default health check and default
   * factory.
   */
  ObjectPool()
      : check_health_([](const DataType &obj) { return true; }),
        factory_([] { return DataType{}; }),
        size_(0) {}

  /**
   * @brief Constructs a pre-sized pool with health check and custom factory.
   * @param size Initial number of objects to pre-create.
   * @param check_health Function to validate objects.
   * @param get_new Factory function to create new DataType instances.
   */
  ObjectPool(size_t size, std::function<bool(const DataType &)> check_health,
             std::function<DataType()> get_new)
      : check_health_(check_health), factory_(get_new), size_(size) {
    for (auto i = 0uz; i < size; i++) {
      auto obj = std::make_unique<DataType>(factory_());
      auto succ = check_health_(*obj);

      while (!succ) {
        obj = std::make_unique<DataType>(factory_());
        succ = check_health_(*obj);
      }

      queue_.push(std::move(obj));
    }
  }

  /**
   * @brief Constructs a pre-sized pool with custom health check and default
   * factory.
   * @param size Initial number of objects to pre-create.
   * @param check_health Function to validate objects.
   */
  ObjectPool(size_t size, std::function<bool(const DataType &)> check_health)
      : check_health_(check_health),
        factory_([] { return DataType{}; }),
        size_(size) {
    for (auto i = 0uz; i < size; i++) {
      auto obj = std::make_unique<DataType>(factory_());
      auto succ = check_health_(*obj);

      while (!succ) {
        obj = std::make_unique<DataType>(factory_());
        succ = check_health_(*obj);
      }

      queue_.push(std::move(obj));
    }
  }

  /**
   * @brief Constructs a pre-sized pool with default health check and custom
   * factory.
   * @param size Initial number of objects to pre-create.
   * @param get_new Factory function to create new DataType instances.
   */
  ObjectPool(size_t size, std::function<DataType()> get_new)
      : check_health_([](const DataType &obj) { return true; }),
        factory_(get_new),
        size_(size) {
    for (auto i = 0uz; i < size; i++) {
      auto obj = std::make_unique<DataType>(factory_());
      auto succ = check_health_(*obj);

      while (!succ) {
        obj = std::make_unique<DataType>(factory_());
        succ = check_health_(*obj);
      }

      queue_.push(std::move(obj));
    }
  }

  /**
   * @brief Constructs a pre-sized pool with default health check and default
   * factory.
   * @param size Initial number of objects to pre-create.
   */
  explicit ObjectPool(size_t size)
      : check_health_([](const DataType &obj) { return true; }),
        factory_([] { return DataType{}; }),
        size_(size) {
    for (auto i = 0uz; i < size; i++) {
      auto obj = std::make_unique<DataType>(factory_());
      auto succ = check_health_(*obj);

      while (!succ) {
        obj = std::make_unique<DataType>(factory_());
        succ = check_health_(*obj);
      }

      queue_.push(std::move(obj));
    }
  }

  ObjectPool(ObjectPool &&) = delete;
  ObjectPool(const ObjectPool &) = delete;
  ObjectPool &operator=(ObjectPool &&) = delete;
  ObjectPool &operator=(const ObjectPool &) = delete;

  /**
   * @class ObjectWrapper
   * @brief RAII wrapper that returns the object to the pool on destruction.
   */
  class ObjectWrapper;

  friend class ObjectWrapper;

  class ObjectWrapper {
    std::unique_ptr<DataType> data_ptr_;   /**< Held object pointer. */
    std::shared_ptr<ObjectPool> pool_ptr_; /**< Owning pool reference. */
    friend class ObjectPool<DataType, Queue>;

    /**
     * @brief Internal constructor used by ObjectPool.
     * @param pool_ptr Shared pointer to the owning pool.
     * @param data_ptr Unique pointer to the pooled object.
     */
    ObjectWrapper(std::shared_ptr<ObjectPool> pool_ptr,
                  std::unique_ptr<DataType> data_ptr)
        : data_ptr_(std::move(data_ptr)), pool_ptr_(pool_ptr_) {}

   public:
    /**
     * @brief Move constructor.
     */
    ObjectWrapper(ObjectWrapper &&old)
        : data_ptr_(std::move(old.data_ptr_)), pool_ptr_(old.pool_ptr_) {}

    // Deleted copy semantics
    ObjectWrapper(const ObjectWrapper &) = delete;

    /**
     * @brief Move assignment operator.
     * @return Reference to this wrapper.
     */
    ObjectWrapper &operator=(ObjectWrapper &&old) {
      this->pool_ptr_ = old.pool_ptr_;
      this->data_ptr_ = std::move(old.data_ptr_);
      return *this;
    }

    // Deleted copy semantics
    ObjectWrapper &operator=(const ObjectWrapper &old) = delete;

    /**
     * @brief Default constructor. Creates an invalid wrapper.
     */
    ObjectWrapper() = default;

    /**
     * @brief Access the managed object pointer.
     * @return Pointer to the DataType instance.
     */
    DataType *operator->() const { return data_ptr_.get(); }

    /**
     * @brief Dereference the managed object.
     * @return Reference to the DataType instance.
     */
    DataType &operator*() const { return *data_ptr_; }

    /**
     * @brief Checks if the wrapper holds a valid object.
     * @return true if valid; false otherwise.
     */
    bool is_valid() {
      return (this->pool_ptr_ != nullptr && this->data_ptr_ != nullptr);
    }

    /**
     * @brief Destructor returns the object to the pool or replaces it if
     * unhealthy.
     */
    ~ObjectWrapper() try {
      bool is_healthy = pool_ptr_->check_health_(*data_ptr_);
      if (is_healthy) {
        pool_ptr_->queue_.push(std::move(data_ptr_));
      } else {
        std::unique_ptr<DataType> obj;
        bool succ = false;

        while (!succ) {
          try {
            obj = std::make_unique<DataType>(pool_ptr_->factory_());
            succ = pool_ptr_->check_health_(*obj);
          } catch (const std::exception &e) {
          }
        }
        pool_ptr_->queue_.push(std::move(obj));
      }
    } catch (const std::exception &e) {
    }
  };

  /**
   * @brief Attempts to pop an object into the provided wrapper.
   * @param[out] obj Wrapper that will hold the object if successful.
   * @return true if an object was obtained and passed health check; false
   * otherwise.
   */
  bool getObject(ObjectWrapper &obj) {
    std::unique_ptr<DataType> obj_ptr;
    auto res = this->queue_.pop(obj_ptr);

    if (res) {
      auto succ = check_health_(*obj_ptr);

      while (!succ) {
        obj_ptr = std::make_unique<DataType>(factory_());
        succ = check_health_(*obj);
      }
    }

    return res;
  }

  /// Only valid when using my safe queue or other compatiable queues.
  /**
   * @brief Pops an object and returns it wrapped, blocking if necessary.
   * @throws std::out_of_range if pool size is zero.
   * @return ObjectWrapper containing a valid object from the pool.
   */
  ObjectWrapper getObject() {
    if (this->size_ == 0) {
      throw std::out_of_range(
          "The size of the pool is zero, please resize first!");
    }

    std::unique_ptr<DataType> obj(this->queue_.pop());
    bool succ = check_health_(*obj);

    if (!succ) {
      obj = std::make_unique<DataType>(factory_());
      succ = check_health_(*obj);
    }

    return ObjectWrapper(this->shared_from_this(), std::move(obj));
  }

  /// The size can only become greator for the case of safety.
  /**
   * @brief Increases the pool size by adding new instances.
   * @param size New desired pool size; must be >= current size.
   */
  void resize(size_t size) {
    auto old_size = this->size_;

    if (size > old_size) {
      auto diff = size - old_size;
      for (auto i = 0uz; i < diff; i++) {
        auto obj = std::make_unique<DataType>(factory_());
        auto succ = check_health_(*obj);

        while (!succ) {
          obj = std::make_unique<DataType>(factory_());
          succ = check_health_(*obj);
        }

        queue_.push(std::move(obj));
      }

      this->size_ = size;
    }
  }

  ~ObjectPool() = default;
};

}  // namespace blib

#endif
