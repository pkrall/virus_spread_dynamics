#pragma once

#include <memory>
#include <mutex>

#include "wstp.h"

namespace services
{
	/* Status code that behaves like an enum class but has methods.
	 * (A real enum class cannot have member functions, so the enum is
	 * wrapped in a class.)
	 *
	 *	MathematicaStatus s = MathematicaProxy::init();
	 *	if (s != MathematicaStatus::OK)
	 *		fprintf(stderr, "%s\n", s.identify());
	 */
	class MathematicaStatus
	{
	public:
		enum Value
		{
			OK = 0,
			WsInitializeFailed,		// WSInitialize failed (or no memory for the proxy)
			LinkOpenFailed,			// WSOpenString failed
			KernelConnectFailed,	// kernel not launched/connected (MathKernel.exe in the path?)
			InvalidInput,			// no valid integer read from stdin
			WstpError,				// WSTP reports an error on the link
			UnexpectedResult		// the kernel's answer does not have the expected form
		};

		constexpr MathematicaStatus(Value value) : value_(value) {}

		// Allows comparisons (==, !=) and switch statements on the status.
		constexpr operator Value() const { return value_; }

		// Prevents accidental use as a condition: if (status) ...
		explicit operator bool() const = delete;

		// Returns the name of the status code, e.g. "LinkOpenFailed".
		// The string is a literal and must not be modified or freed.
		const char* identify() const;

	private:
		Value value_;
	};

	/* Singleton that launches the Mathematica kernel and keeps the WSTP link.
	 *	MathKernel.exe must be in the path!
	 *
	 * Usage:
	 *	if (services::MathematicaProxy::init() == services::MathematicaStatus::OK)
	 *		services::MathematicaProxy::instance()->test_mathematica(&n);
	 */
	class MathematicaProxy
	{
	public:
		// Creates the singleton (launches the kernel). If it already exists,
		// nothing happens and OK is returned. If construction fails, the error
		// code is returned and instance() stays nullptr; init() can then be
		// called again.
		static MathematicaStatus init();

		// Returns the singleton, or nullptr if init() has not been called
		// successfully yet.
		static MathematicaProxy* instance();

		// Factors *integer_to_test and prints the prime factors.
		// If nullptr is passed, the integer is read from stdin.
		MathematicaStatus test_mathematica(int* integer_to_test);

		// Closes the link and thereby terminates the kernel.
		// Public because std::unique_ptr has to call the destructor.
		~MathematicaProxy();

		MathematicaProxy(const MathematicaProxy&) = delete;
		MathematicaProxy& operator=(const MathematicaProxy&) = delete;

	private:
		MathematicaProxy();

		// Reports the error on stderr, resets the link and returns the
		// matching error code.
		MathematicaStatus fail();

		static std::unique_ptr<MathematicaProxy> instance_;
		static std::mutex init_mutex_;

		MathematicaStatus status_ = MathematicaStatus::OK;	// result of construction
		WSENV ep_ = (WSENV)0;
		WSLINK lp_ = (WSLINK)0;
		std::mutex mutex_;	// a WSTP link is not thread-safe
	};
}
