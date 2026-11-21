#include "pch.h"
#include "App.xaml.h"
#include <winrt/Microsoft.UI.Xaml.h>
#include "Services/Mathematica/MathematicaProxy.h"
#include <iostream>
#include <conio.h>

using namespace winrt::TRTU_TRTA_CSCP::implementation;


int main(int argc, char* argv[])
{
    // first check whether we have required services
    services::MathematicaStatus mathematica_status = services::MathematicaStatus::OK;
    if (   (mathematica_status = services::MathematicaProxy::init()) != services::MathematicaStatus::OK
        || (mathematica_status = services::MathematicaProxy::instance()->test_mathematica(nullptr)) != services::MathematicaStatus::OK) {
        printf("failed to and test start mathematica, status=%s", mathematica_status.identify());
        char c = _getch();
        exit(-1);
    }

    winrt::init_apartment(winrt::apartment_type::single_threaded);

    winrt::Microsoft::UI::Xaml::Application::Start([](auto&&)
        {
            // Hier wird die App-Instanz erzeugt und das MainWindow gestartet
            winrt::make<App>();
        });


    return 0;
}