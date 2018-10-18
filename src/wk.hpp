/*-
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Guram Duka
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
//------------------------------------------------------------------------------
#ifndef WK_HPP_INCLUDED
#define WK_HPP_INCLUDED
//------------------------------------------------------------------------------
#include <functional>
#include <thread>
#include <queue>
#include <mutex>
#include <shared_mutex>
#include <memory>
#include <condition_variable>
//------------------------------------------------------------------------------
namespace nn {  // namespace NaturalNumbers
//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
struct cdc256 {
	struct data {
		uint64_t a, b, c, d, e, f, g, h;

		void shuffle() {
			a -= e; f ^= h >>  9; h += a;
			b -= f; g ^= a <<  9; a += b;
			c -= g; h ^= b >> 23; b += c;
			d -= h; a ^= c << 15; c += d;
			e -= a; b ^= d >> 14; d += e;
			f -= b; c ^= e << 20; e += f;
			g -= c; d ^= f >> 17; f += g;
			h -= d; e ^= g << 14; g += h;
		}

		void update(const data & v) {
			a -= v.e; f ^= v.h >>  9; h += v.a;
			b -= v.f; g ^= v.a <<  9; a += v.b;
			c -= v.g; h ^= v.b >> 23; b += v.c;
			d -= v.h; a ^= v.c << 15; c += v.d;
			e -= v.a; b ^= v.d >> 14; d += v.e;
			f -= v.b; c ^= v.e << 20; e += v.f;
			g -= v.c; d ^= v.f >> 17; f += v.g;
			h -= v.d; e ^= v.g << 14; g += v.h;
		}
	};

	data data_;
	uint64_t position_;

	void init() {
		static const data init_value = {
			UINT64_C(0xA640524A5B44F1FC),
			UINT64_C(0xC535059705F0BB7E),
			UINT64_C(0xC8ED76CF6B6EA626),
			UINT64_C(0x531D1E8E254EA59E),
			UINT64_C(0x8C0FE7F3E46E2A80),
			UINT64_C(0x1C53F41FD1E3A7F8),
			UINT64_C(0x08D4DEAAA1C33335),
			UINT64_C(0x4C592980FBE9B011)
		};
		data_ = init_value;
		position_ = 0;
	}

	void update(const void * p, uintptr_t size) {
		data pad;

		auto dig = [] (data & d, const data & p) {
			d.update(p);
			d.shuffle();
		};

		position_ += size;

		while( size >= sizeof(data) ) {
			dig(data_, *reinterpret_cast<const data *>(p));
			p = (const uint8_t *) p + sizeof(data);
			size -= sizeof(data);
		}

		if( size > 0 ) {
			memcpy(&pad, p, size);
			memset((uint8_t *) &pad + size, 0, sizeof(data) - size);
			dig(data_, pad);
		}

		pad.a = pad.b = pad.c = pad.d = pad.e = pad.f = pad.g = pad.h = position_;
		dig(data_, pad);
	}

	template <typename T>
	void update(T & container) {
		update(container.data(), container.size());
	}

	//std::array<uint8_t, 64> final() {
	//	return *reinterpret_cast<std::array<uint8_t, 64> *>(&data_);
	//}
};
//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
struct cdc512 {
	struct data {
		uint64_t a, b, c, d, e, f, g, h;
		uint64_t i, j, k, l, m, n, o, p;

		void shuffle() {
			a -= i; i ^= h >> 10; p += a;
			b -= j; j ^= g << 10; o += b;
			c -= k; k ^= f >> 23; n += c;
			d -= l; l ^= e << 15; m += d;
			e -= m; m ^= d >> 14; l += e;
			f -= n; n ^= c << 20; k += f;
			g -= o; o ^= b >> 17; j += g;
			h -= p; p ^= a << 14; i += h;
			i -= a; a ^= p >> 10; h += i;
			j -= b; b ^= o << 10; g += j;
			k -= c; c ^= n >> 23; f += k;
			l -= d; d ^= m << 15; e += l;
			m -= e; e ^= l >> 14; d += m;
			n -= f; f ^= k << 20; c += n;
			o -= g; g ^= j >> 17; b += o;
			p -= h; h ^= i << 14; a += p;
		}

		void update(const data & v) {
			a -= v.i; i ^= v.h >> 10; p += v.a;
			b -= v.j; j ^= v.g << 10; o += v.b;
			c -= v.k; k ^= v.f >> 23; n += v.c;
			d -= v.l; l ^= v.e << 15; m += v.d;
			e -= v.m; m ^= v.d >> 14; l += v.e;
			f -= v.n; n ^= v.c << 20; k += v.f;
			g -= v.o; o ^= v.b >> 17; j += v.g;
			h -= v.p; p ^= v.a << 14; i += v.h;
			i -= v.a; a ^= v.p >> 10; h += v.i;
			j -= v.b; b ^= v.o << 10; g += v.j;
			k -= v.c; c ^= v.n >> 23; f += v.k;
			l -= v.d; d ^= v.m << 15; e += v.l;
			m -= v.e; e ^= v.l >> 14; d += v.m;
			n -= v.f; f ^= v.k << 20; c += v.n;
			o -= v.g; g ^= v.j >> 17; b += v.o;
			p -= v.h; h ^= v.i << 14; a += v.p;
		}
	};

	data data_;
	uint64_t position_;

	void init() {
		static const data init_value = {
			UINT64_C(0xA640524A5B44F1FC),
			UINT64_C(0xC535059705F0BB7E),
			UINT64_C(0xC8ED76CF6B6EA626),
			UINT64_C(0x531D1E8E254EA59E),
			UINT64_C(0x8C0FE7F3E46E2A80),
			UINT64_C(0x1C53F41FD1E3A7F8),
			UINT64_C(0x08D4DEAAA1C33335),
			UINT64_C(0x4C592980FBE9B011),
			UINT64_C(0x992E367BE6F0EA1E),
			UINT64_C(0x71DCF41FFACC283F),
			UINT64_C(0xC9581F48D85ABD75),
			UINT64_C(0xE4B93335FF1CE990),
			UINT64_C(0xE51D6424EFEC1E01),
			UINT64_C(0x353867A0E66C2A39),
			UINT64_C(0xA8DBF7B782226B67),
			UINT64_C(0x9F8B7F0DC254488E)
		};
		data_ = init_value;
		position_ = 0;
	}

	void update(const void * p, uintptr_t size) {
		data pad;

		auto dig = [] (data & d, const data & p) {
			d.update(p);
			d.shuffle();
		};

		position_ += size;

		while( size >= sizeof(data) ) {
			dig(data_, *reinterpret_cast<const data *>(p));
			p = (const uint8_t *) p + sizeof(data);
			size -= sizeof(data);
		}

		if( size > 0 ) {
			memcpy(&pad, p, size);
			memset((uint8_t *) &pad + size, 0, sizeof(data) - size);
			dig(data_, pad);
		}

		pad.a = pad.b = pad.c = pad.d = pad.e = pad.f = pad.g = pad.h = position_;
		dig(data_, pad);
	}

	template <typename T>
	void update(T & container) {
		update(container.data(), container.size());
	}

	//std::array<uint8_t, 128> final() {
	//	return *reinterpret_cast<std::array<uint8_t, 128> *>(&data_);
	//}
};
//------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
template <typename FnType>
class ThreadPoolT {
	public:
		typedef std::function<FnType> TaskType;

		~ThreadPoolT() {
			shutdown_ = true;

			for( auto & a : cvs_ )
				a->notify_all();

			for( auto & i : threads_ )
				i->join();
		}

		ThreadPoolT(uintptr_t threads = 0) :
		    threads_(),
		    cvs_(),
		    runs_(),
		    wv_(),
		    mutex_(),
		    tasks_(),
		    counter_(0),
		    done_(false),
		    shutdown_(false) {

			if( threads == 0 )
				threads = std::thread::hardware_concurrency();

			tasks_.reserve(threads * 2);

			for( uintptr_t i = 0; i < threads; i++ ) {
				cvs_.emplace_back(new std::condition_variable);
				runs_.emplace_back(false);
			}

			while( threads-- > 0 ) {
				threads_.emplace_back(new std::thread(&ThreadPoolT<FnType>::run, this, threads_.size()));
			}
		}

		void run_task(const TaskType & task) {
			//std::unique_lock<decltype(mutex_)> locker(mutex_);
			tasks_.emplace_back(task);
		}

		void wait() {
			counter_ = tasks_.size() >= threads_.size() ? threads_.size() : tasks_.size();

			if( counter_ ) {
				std::unique_lock<decltype(mutex_)> locker(mutex_);
				for( uintptr_t i = 0; i < counter_; i++ ) {
					runs_[i] = true;
					cvs_[i]->notify_one();
				}
				wv_.wait(locker, [&] { return done_; });
				tasks_.clear();
			}

			done_ = false;
		}

	private:
		std::vector<std::unique_ptr<std::thread>> threads_;
		std::vector<std::unique_ptr<std::condition_variable>> cvs_;
		std::vector<bool> runs_;
		std::condition_variable wv_;
		std::mutex mutex_;
		std::vector<TaskType> tasks_;
		uintptr_t counter_;
		bool done_;
		bool shutdown_;

		void run(uintptr_t index) {
			std::unique_lock<decltype(mutex_)> locker(mutex_);

			while( !shutdown_ ) {
				cvs_[index]->wait(locker, [&] { return runs_[index] || shutdown_; });

				if( runs_[index] ) {
					locker.unlock();

					uintptr_t a = 0;

					for( auto i = index; i < tasks_.size(); i += threads_.size() ) {
						tasks_[i]();
						a = 1;
					}

					locker.lock();

					if( a && --counter_ == 0 ) {
						done_ = true;
						wv_.notify_all();
					}

					runs_[index] = false;
				}
			}
		}
};
//------------------------------------------------------------------------------
typedef ThreadPoolT<void()> ThreadPool;
//------------------------------------------------------------------------------
} // namespace NaturalNumbers
//------------------------------------------------------------------------------
#endif // WK_HPP_INCLUDED
//------------------------------------------------------------------------------
