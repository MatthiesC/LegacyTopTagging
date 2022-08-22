#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include <iostream>

using namespace std;
using namespace uhh2;

class DummyAnalysisModule: public AnalysisModule {
public:
    DummyAnalysisModule(Context &){}
    virtual bool process(Event & event){
        cout << "Event: " << event_i++ << endl;
        int n_ele = 0;
        for(const Electron & ele : *event.electrons) {
            cout << "Electron " << n_ele++ << endl;
            cout << to_string(ele.energy()) << endl;
            cout << to_string(ele.get_tag(Electron::tag::residual_ecalTrkEnergyPostCorr)) << endl;
            cout << to_string(ele.get_tag(Electron::tag::residual_energyScaleUp)) << endl;
            cout << to_string(ele.get_tag(Electron::tag::residual_energyScaleDown)) << endl;
        }
        int n_pho = 0;
        for(const Photon & pho : *event.photons) {
            cout << "Photon " << n_pho++ << endl;
            cout << to_string(pho.energy()) << endl;
            cout << to_string(pho.get_tag(Photon::tag::residual_ecalEnergyPostCorr)) << endl;
            cout << to_string(pho.get_tag(Photon::tag::residual_energyScaleUp)) << endl;
            cout << to_string(pho.get_tag(Photon::tag::residual_energyScaleDown)) << endl;
        }
        return true;
    }
    virtual void endInputData(){}
    long event_i = 0;
};

UHH2_REGISTER_ANALYSIS_MODULE(DummyAnalysisModule)
