#include "pch.h"
#include <stdio.h>
#include <new>

#include "MathematicaProxy.h"

namespace services
{
	const char* MathematicaStatus::identify() const
	{
		switch (value_)
		{
		case OK:					return "OK";
		case WsInitializeFailed:	return "WsInitializeFailed";
		case LinkOpenFailed:		return "LinkOpenFailed";
		case KernelConnectFailed:	return "KernelConnectFailed";
		case InvalidInput:			return "InvalidInput";
		case WstpError:				return "WstpError";
		case UnexpectedResult:		return "UnexpectedResult";
		}
		return "Unknown";
	}


	std::unique_ptr<MathematicaProxy> MathematicaProxy::instance_;
	std::mutex MathematicaProxy::init_mutex_;


	MathematicaStatus MathematicaProxy::init()
	{
		std::lock_guard<std::mutex> lock(init_mutex_);

		if (instance_)
			return MathematicaStatus::OK;

		// No make_unique, because the constructor is private.
		instance_.reset(new (std::nothrow) MathematicaProxy());
		if (!instance_)
			return MathematicaStatus::WsInitializeFailed;

		MathematicaStatus status = instance_->status_;
		if (status != MathematicaStatus::OK)
			instance_.reset();	// the destructor cleans up whatever was already opened

		return status;
	}


	MathematicaProxy* MathematicaProxy::instance()
	{
		return instance_.get();
	}


	MathematicaProxy::MathematicaProxy()
	{
		int err = 0;

		ep_ = WSInitialize((WSParametersPointer)0);
		if (ep_ == (WSENV)0) {
			status_ = MathematicaStatus::WsInitializeFailed;
			return;
		}

		lp_ = WSOpenString(ep_, " -linkmode launch -linkname MathKernel.exe", &err);
		if (lp_ == (WSLINK)0) {
			fprintf(stderr, "WSOpenString failed (WSTP error %d).\n", err);
			status_ = MathematicaStatus::LinkOpenFailed;
			return;
		}

		// Make sure the kernel has really been launched and is connected.
		if (!WSActivate(lp_)) {
			fprintf(stderr, "Could not connect to MathKernel: %s.\n", WSErrorMessage(lp_));
			status_ = MathematicaStatus::KernelConnectFailed;
			return;
		}
	}


	MathematicaProxy::~MathematicaProxy()
	{
		// Closing the link terminates the kernel that was started in "launch" mode.
		if (lp_) WSClose(lp_);
		if (ep_) WSDeinitialize(ep_);
	}


	MathematicaStatus MathematicaProxy::test_mathematica(int* integer_to_test)
	{
		int pkt, n, prime, expt;
		int len, lenp, k;

		if (integer_to_test == nullptr)
		{
			printf("Integer to factor: ");

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
			if (scanf_s("%d", &n) != 1)
#else
			if (scanf("%d", &n) != 1)
#endif
				return MathematicaStatus::InvalidInput;
		}
		else {
			n = *integer_to_test;
		}

		std::lock_guard<std::mutex> lock(mutex_);

		WSPutFunction(lp_, "EvaluatePacket", 1L);
		WSPutFunction(lp_, "FactorInteger", 1L);
		WSPutInteger(lp_, n);
		WSEndPacket(lp_);

		while ((pkt = WSNextPacket(lp_), pkt) && pkt != RETURNPKT) {
			WSNewPacket(lp_);
			if (WSError(lp_)) return fail();
		}

		if (!WSTestHead(lp_, "List", &len)) return fail();
		for (k = 1; k <= len; k++) {
			if (WSTestHead(lp_, "List", &lenp)
				&& lenp == 2
				&& WSGetInteger(lp_, &prime)
				&& WSGetInteger(lp_, &expt)
				) {
				printf("%d ^ %d\n", prime, expt);
			}
			else {
				return fail();
			}
		}

		// Finish the packet so the link is ready for the next call.
		WSNewPacket(lp_);

		return MathematicaStatus::OK;
	}


	MathematicaStatus MathematicaProxy::fail()
	{
		MathematicaStatus status = MathematicaStatus::UnexpectedResult;

		if (WSError(lp_)) {
			fprintf(stderr, "Error detected by WSTP: %s.\n", WSErrorMessage(lp_));
			status = MathematicaStatus::WstpError;
		}
		else {
			fprintf(stderr, "Error detected by this program.\n");
		}

		// Try to keep the link usable for further calls.
		WSClearError(lp_);
		WSNewPacket(lp_);

		return status;
	}
}
