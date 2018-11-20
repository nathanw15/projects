#define GAMMA_H_INC_ALL         // define this to include all header files
#define GAMMA_H_NO_IO           // define this to avoid bringing AudioIO from Gamma

#include "Gamma/Gamma.h"
#include "al/core.hpp"
#include "al/core/sound/al_Speaker.hpp"
#include "al/util/ui/al_Parameter.hpp"
#include "al/util/ui/al_ParameterServer.hpp"
#include "Gamma/SamplePlayer.h"
#include "al/util/ui/al_ControlGUI.hpp"

#include <atomic>
#include <vector>

#define SAMPLE_RATE 44100
#define BLOCK_SIZE (2048)

using namespace al;
using namespace std;

ParameterServer paramServer("127.0.0.1",8080);
Parameter srcAzimuth("srcAzimuth","",-0.5,"",-1.0*M_PI,M_PI);
ParameterBool updatePanner("updatePanner","",0.0);

//ParameterBool useRamp("useRamp","",0);
//Parameter rampStartAzimuth("rampStartAzimuth","",-0.5,"",-1.0*M_PI,M_PI);
//Parameter rampEndAzimuth("rampEndAzimuth","",0.5,"",-1.0*M_PI,M_PI);
//Parameter rampDuration("rampDuration", "",1.0,"",0.0,5.0);

float radius = 5.0;

Parameter sourceSound("sourceSound","",0.0,"",0.0,3.0);
Parameter soundFileIdx("soundFileIdx","",0.0,"",0.0,3.0);

ParameterBool sampleWise("sampleWise","",0.0);
ParameterBool useDelay("useDelay","", 0.0);

Parameter masterGain("masterGain","",0.5,"",0.0,1.0);
ParameterBool soundOn("soundOn","",0.0);

//ParameterBool triggerRamp("triggerRamp","",0.0);

//int sourceSound = 1;
//int soundFileIdx = 1;

class VirtualSource {
public:

    float aziInRad;
    gam::SamplePlayer<> samplePlayer;

    int sound = 0;
    int fileIdx = 0;

    VirtualSource(){

    }

    void onSoundSample(AudioIOData &io){

    }

};

class SpeakerV: public Speaker {
public:

    ParameterBool *enabled;
    std::string oscTag;
    float aziInRad;

    int delay;
    void *buffer;
    int bufferSize;
    int readPos,writePos;

    SpeakerV(int chan, float az=0.f, float el=0.f, int gr=0, float rad=1.f, float ga=1.f, int del = 0){
        delay = del;
        readPos = 0;
        writePos = delay;
//        bufferSize = BLOCK_SIZE+(delay*2);
        bufferSize = 44100*2;

        buffer = calloc(bufferSize,sizeof(float));
        //cout << "Speaker: " << chan << " Write Pos: " << writePos << endl;
         //buffer = new float[bufferSize]();

      // buffer = new gam::Array<float>();
      // buffer.resize(delay*2);

        deviceChannel = chan;
        azimuth= az;
        elevation = el;
        group = gr;
        radius = rad;
        gain = ga;
        aziInRad = toRad(az);

        oscTag = "speaker"+ std::to_string(deviceChannel) + "/enabled";
        enabled = new ParameterBool(oscTag,"",1.0);
        enabled->setHint("latch",1.f);
        paramServer.registerParameter(*enabled);
    }

    float read(){
        if(readPos >= bufferSize){
            readPos = 0;
        }
        float *b = (float*)buffer;
        float val = b[readPos];
        readPos++;
        return val;
    }

    void write(float samp){
        if(writePos >= bufferSize){
            writePos = 0;
        }
        float *b = (float*)buffer;
        b[writePos] = samp;
        writePos++;
    }
};

std::vector<SpeakerV> speakers;
std::vector<SpeakerV*> enabledSpeakers;
Mat<2,double> matrix;


struct Ramp {

    unsigned int startSample, endSample;
    bool done = false;
    bool running = false;

    //Parameter *trigger;//("triggerRamp","",0.0);
    bool trigger = false;

    ParameterBundle rampBundle{"rampBundle"};

    ParameterBool useRamp{"useRamp","",0};
    ParameterBool triggerRamp{"triggerRamp","",0.0};
    Parameter rampStartAzimuth{"rampStartAzimuth","",-0.5,"",-1.0*M_PI,M_PI};
    Parameter rampEndAzimuth{"rampEndAzimuth","",0.5,"",-1.0*M_PI,M_PI};
    Parameter rampDuration{"rampDuration", "",1.0,"",0.0,5.0};

    Ramp(){
        //trigger = new Parameter("triggerRamp","",0.0);
        //paramServer.registerParameter(*trigger);
        useRamp.setHint("latch",1.f);
        rampBundle << useRamp << triggerRamp << rampStartAzimuth << rampEndAzimuth << rampDuration;

        triggerRamp.registerChangeCallback([&](float val){
            if(val == 1.f){ // is this correct way to check?
                //cout << "Ramp triggered" << endl;
                trigger = true;
//                gam::SamplePlayer<> *player = samplePlayers[soundFileIdx];
//                player->reset();
            }
        });
    }

    void set(float startAzi, float endAzi, float dur){
//        startAzimuth = startAzi;
//        endAzimuth = endAzi;
//        duration = dur;
        rampStartAzimuth.set(startAzi);
        rampEndAzimuth.set(endAzi);
        rampDuration.set(dur);
    }

    void start(unsigned int startSamp){
        //unsigned int i=static_cast<unsigned int>(duration);
        startSample = startSamp;
        endSample = startSamp +  rampDuration.get() *SAMPLE_RATE;
        running = true;
        done = false;
        //cout << startSamp << " endSamp: " << endSample << endl;
    }

    float next(unsigned int sampleNum){

        if(trigger){
            trigger = false;
            start(sampleNum);
        }

        if(!done && !running){
            return rampStartAzimuth.get();
        }else if(!done && running){

            if(sampleNum > endSample){
                sampleNum = endSample;
                done = true;
                running = false;
            }
            float val = (((sampleNum - startSample) * (rampEndAzimuth.get() - rampStartAzimuth.get())) / (endSample - startSample)) + rampStartAzimuth.get();
            return val;
        } else {
            return rampEndAzimuth.get();
        }
    }
};

//bool speakerSort(SpeakerV const &first, SpeakerV const &second){
//    return first.azimuth < second.azimuth;
//}

bool speakerSort(SpeakerV const *first, SpeakerV const *second){
    return first->azimuth < second->azimuth;
}

void initPanner(){
    enabledSpeakers.clear();
    for(int i = 0; i < speakers.size(); i ++){
        if(speakers[i].enabled->get() > 0.5){
            enabledSpeakers.push_back(&speakers[i]);
        }
    }
    std::sort(enabledSpeakers.begin(),enabledSpeakers.end(),&speakerSort);
}

class MyApp : public App
{
public:
    Mesh mSpeakerMesh;
    vector<Mesh> mVec;
    vector<int>  sChannels;
    SpeakerLayout speakerLayout;
    atomic<float> *mPeaks {nullptr};
    float speedMult = 0.03f;
    Vec3d srcpos {0.0,0.0,0.0};
    Ramp linearRamp;
    vector<gam::SamplePlayer<>*> samplePlayers;
    SearchPaths searchpaths;
    ControlGUI parameterGUI;

    vector<VirtualSource> sources;

    MyApp()
    {
//        parameterGUI << soundOn << srcAzimuth << updatePanner << triggerRamp << rampStartAzimuth << rampEndAzimuth << rampDuration << sourceSound << soundFileIdx << sampleWise << useDelay << masterGain << useRamp;
        parameterGUI << soundOn << srcAzimuth << updatePanner << linearRamp.rampBundle << sourceSound << soundFileIdx << sampleWise << useDelay << masterGain;


        //        parameterServer() << srcAzimuth << updatePanner << rampStartAzimuth << rampEndAzimuth << rampDuration << sourceSound << soundFileIdx << sampleWise;

       // paramServer << srcAzimuth << updatePanner << triggerRamp << rampStartAzimuth << rampEndAzimuth << rampDuration << sourceSound << soundFileIdx << sampleWise << useDelay << masterGain << useRamp << soundOn;
        soundOn.setHint("latch", 1.f);
        sampleWise.setHint("latch", 1.f);
        useDelay.setHint("latch",1.f);

        sourceSound.setHint("intcombo",1.f);
        soundFileIdx.setHint("intcombo",1.f);

        updatePanner.registerChangeCallback([&](float val){
            if(val == 1.f){ // is this correct way to check?
                cout << "panner Updated" << endl;
            }
            initPanner();
        });
    }

    void createMatrix(Vec3d left, Vec3d right){
        matrix.set(left.x,left.y,right.x,right.y);
    }

    Vec3d ambiSphericalToOGLCart(float azimuth, float radius){
        Vec3d ambiSpherical;
        float elevation = 0.0;

        //find xyz in cart audio coords
        float x = radius * cos(elevation) * cos(azimuth);
        float y = radius * cos(elevation) * sin(azimuth);
        float z = radius * sin(elevation);

        //convert to open_gl cart coords
        ambiSpherical[0] = -y;
        ambiSpherical[1] = z;
        ambiSpherical[2] = -x;

        return ambiSpherical;
    }

    void openGLCartToAmbiCart(Vec3f &vec){
        Vec3f tempVec = vec;
        vec[0] = tempVec.z*-1.0;
        vec[1] = tempVec.x*-1.0;
        vec[2] = tempVec.y;
    }

    //TODO: io not used here
    Vec3d calcGains(AudioIOData &io, const float &srcAzi, int &speakerChan1, int &speakerChan2){

        srcpos = ambiSphericalToOGLCart(srcAzi,radius);
        Vec3f ambiCartSrcPos = srcpos;
        openGLCartToAmbiCart(ambiCartSrcPos);
        std::sort(enabledSpeakers.begin(),enabledSpeakers.end(),&speakerSort);
        Vec3d gains(0.,0.,0.);
        float speakSrcAngle,linearDistance;

        //check if source is beyond the first or last speaker
        if(srcAzi < enabledSpeakers[0]->aziInRad){
            speakerChan1 = enabledSpeakers[0]->deviceChannel;
            speakerChan2 = enabledSpeakers[0+1]->deviceChannel;
            speakSrcAngle = fabsf(enabledSpeakers[0]->aziInRad - srcAzi);
            gains.x = 1.f / radius * (M_PI - speakSrcAngle);

        } else if(srcAzi > enabledSpeakers[enabledSpeakers.size()-1]->aziInRad){
            speakerChan1 = enabledSpeakers[enabledSpeakers.size()-2]->deviceChannel;//set to speaker before last
            speakerChan2 = enabledSpeakers[enabledSpeakers.size()-1]->deviceChannel;
            speakSrcAngle = fabsf(enabledSpeakers[enabledSpeakers.size()-1]->aziInRad - srcAzi);
            linearDistance = 2.0*radius*cos((M_PI - speakSrcAngle)/2.0);
            gains.y = 1.f / radius * (M_PI - speakSrcAngle);

        } else{//Source between first and last speakers
            for(int i = 0; i < enabledSpeakers.size()-1; i++){
                speakerChan1 = enabledSpeakers[i]->deviceChannel;
                speakerChan2 = enabledSpeakers[i+1]->deviceChannel;
                if(srcAzi == enabledSpeakers[i]->aziInRad ){
                    gains.x = 1.0;
                    break;
                }else if(srcAzi > enabledSpeakers[i]->aziInRad && srcAzi < enabledSpeakers[i+1]->aziInRad){
                    createMatrix(enabledSpeakers[i]->vec(),enabledSpeakers[i+1]->vec());
                    invert(matrix);
                    for (unsigned i = 0; i < 2; i++){
                        for (unsigned j = 0; j < 2; j++){
                            gains[i] += ambiCartSrcPos[j] * matrix(j,i);
                        }
                    }
                    gains = gains.normalize();
                    break;
                } else if(srcAzi == enabledSpeakers[i+1]->aziInRad){
                    gains.y = 1.0;
                    break;
                }
            }
        }
        return gains;
    }

    void renderBufferDelaySpeakers(AudioIOData &io,const float &srcAzi, const float *buffer){
        int speakerChan1, speakerChan2;
        Vec3d gains = calcGains(io,srcAzi, speakerChan1, speakerChan2);

        for(int i = 0; i < enabledSpeakers.size(); i++){
            SpeakerV *s = enabledSpeakers[i];
            for(int j = 0; j < io.framesPerBuffer();j++){
                if(s->deviceChannel == speakerChan1){
                    s->write(buffer[j]*gains[0]);
                }else if(s->deviceChannel == speakerChan2){
                    s->write(buffer[j]*gains[1]);
                }else{
                    s->write(0.0);
                }
                io.out(s->deviceChannel,j) = s->read();
            }
        }
    }

    void renderSampleDelaySpeakers(AudioIOData &io,const float &srcAzi, const float &sample){
        int speakerChan1, speakerChan2;
        Vec3d gains = calcGains(io,srcAzi, speakerChan1, speakerChan2);

        for(int i = 0; i < enabledSpeakers.size(); i++){
            SpeakerV *s = enabledSpeakers[i];

            if(s->deviceChannel == speakerChan1){
                s->write(sample*gains[0]);
            }else if(s->deviceChannel == speakerChan2){
                s->write(sample*gains[1]);
            }else{
                s->write(0.0);
            }
            io.out(s->deviceChannel,io.frame()) = s->read();
        }
    }


    void renderBuffer(AudioIOData &io,const float &srcAzi, const float *buffer){
        int speakerChan1, speakerChan2;
        Vec3d gains = calcGains(io,srcAzi, speakerChan1, speakerChan2);
        for(int i = 0; i < io.framesPerBuffer(); i++){
            io.out(speakerChan1,i) += buffer[i]*gains[0];
            io.out(speakerChan2,i) += buffer[i]*gains[1];
        }
    }

    void renderSample(AudioIOData &io, const float &srcAzi, const float &sample){
        int speakerChan1, speakerChan2;
        Vec3d gains = calcGains(io,srcAzi, speakerChan1, speakerChan2);
        io.out(speakerChan1,io.frame()) += sample * gains[0];
        io.out(speakerChan2,io.frame()) += sample * gains[1];
    }

    void onInit() override {

        searchpaths.addAppPaths();
        searchpaths.addRelativePath("../sounds");

        samplePlayers.push_back(new gam::SamplePlayer<>(searchpaths.find("lowBoys.wav").filepath().c_str()));
        samplePlayers.push_back(new gam::SamplePlayer<>(searchpaths.find("devSamples.wav").filepath().c_str()));

        linearRamp.set(-2.0,2.0,linearRamp.rampDuration.get());


       // SpeakerV(0,90.0,0.0,0,5.0,1.0,0);
       // speakers.push_back(SpeakerV(0,90.0 - (10.0*0),0.0,0,5.0,1.0,0));

        for (int i = 0; i < 20; i++){
            int delay = rand() % static_cast<int>(44100 + 1);
            speakers.push_back(SpeakerV(i,90.0 - (10.0*i),0.0,0,5.0,0,delay));
        }



//speakers.push_back(SpeakerV(0,90.0,0.0,0,5.0,0,0));
//speakers.push_back(SpeakerV(1,0.0,0.0,0,5.0,0,5000));

//        speakers.push_back(SpeakerV(0,-45.0,0.0,0,5.0));
//        speakers.push_back(SpeakerV(1,-20.0,0.0,0,5.0));
//        speakers.push_back(SpeakerV(2,0.0,0.0,0,5.0));
//        speakers.push_back(SpeakerV(3,20.0,0.0,0,5.0));
//        speakers.push_back(SpeakerV(4,45.0,0.0,0,5.0));

        //speakers[1].enabled->set(0.0);

        initPanner();

        for(int i = 0; i < speakers.size(); i++){
            parameterGUI << speakers[i].enabled;
            //speakers[i].enabled->setHint("latch",1.f);
        }

        int highestChannel = 0;
        for(SpeakerV s:speakers){
            if((int) s.deviceChannel > highestChannel){
                highestChannel = s.deviceChannel;
            }
        }

        audioIO().close();
        audioIO().channelsOut(highestChannel + 1);
        audioIO().open();

        mPeaks = new atomic<float>[highestChannel + 1];

        addSphere(mSpeakerMesh, 1.0, 5, 5);

    }

    void onCreate() override {
        nav().pos(0, 1, 20);
        parameterGUI.init();
    }

    void onAnimate( double dt) override {
        navControl().active(!parameterGUI.usingInput());
    }

    virtual void onSound(AudioIOData &io) override {

        static unsigned int t = 0;
        double sec;
        float srcBuffer[BLOCK_SIZE];

        float mGain = masterGain.get();
        gam::Sync::master().spu(audioIO().fps());

        if(soundOn.get() > 0.5){
            while (io()) {
                int i = io.frame();
                if(sourceSound.get() == 0){
                    float env = (22050 - (t % 22050))/22050.0;
                    srcBuffer[i] = mGain * rnd::uniform() * env;
                } else if(sourceSound.get() ==1){
                    gam::SamplePlayer<> *player = samplePlayers[soundFileIdx];
                    if(player->done()){
                        player->reset();
                    }
                    srcBuffer[i] = mGain * player->operator ()();
                }
                ++t;

                if(linearRamp.useRamp.get()){
                    srcAzimuth.set(linearRamp.next(t));
                }

                if(sampleWise.get() == 1.f){
                    if(useDelay.get() == 1.f){
                        renderSampleDelaySpeakers(io,srcAzimuth.get(),srcBuffer[i]);
                    } else {
                        renderSample(io,srcAzimuth.get(), srcBuffer[i]);
                    }
                }
            }

            if(sampleWise.get() == 0.f){
                if(useDelay.get() == 1.f){
                    renderBufferDelaySpeakers(io,srcAzimuth.get(), srcBuffer);
                }else{
                    renderBuffer(io,srcAzimuth.get(), srcBuffer);
                }
            }
        }

        for (int i = 0; i < speakers.size(); i++) {
            mPeaks[i].store(0);
        }
        for (int speaker = 0; speaker < speakers.size(); speaker++) {
            float rms = 0;
            int deviceChannel = speakers[speaker].deviceChannel;
            for (int i = 0; i < io.framesPerBuffer(); i++) {
                if(deviceChannel < io.channelsOut()) {
                    float sample = io.out(deviceChannel, i);
                    rms += sample * sample;
                }
            }
            rms = sqrt(rms/io.framesPerBuffer());
            mPeaks[deviceChannel].store(rms);
        }
    }

    virtual void onDraw(Graphics &g) override {

        g.clear(0);
        g.blending(true);
        g.blendModeAdd();

        //Draw the source
        g.pushMatrix();
        g.translate(srcpos);
        g.scale(0.3);
        g.color(0.4,0.4, 0.4, 0.5);
        g.draw(mSpeakerMesh);
        g.popMatrix();

        // Draw line
        Mesh lineMesh;
        lineMesh.vertex(0.0,0.0, 0.0);
        lineMesh.vertex(srcpos.x,0.0, srcpos.z);
        lineMesh.vertex(srcpos);
        lineMesh.index(0);
        lineMesh.index(1);
        lineMesh.index(1);
        lineMesh.index(2);
        lineMesh.index(2);
        lineMesh.index(0);
        lineMesh.primitive(Mesh::LINES);
        g.color(1);
        g.draw(lineMesh);

        //Draw the speakers
        for(int i = 0; i < speakers.size(); ++i){
            int devChan = speakers[i].deviceChannel;
            g.pushMatrix();
            g.translate(speakers[i].vecGraphics());
            float peak = mPeaks[devChan].load();
            g.scale(0.02 + peak * 6);
            g.polygonLine();

            if(devChan == 0){
                g.color(1.0,0.0,0.0);
            }else{
            g.color(1);
            }
            g.draw(mSpeakerMesh);
            g.popMatrix();
        }

        parameterGUI.draw(g);
    }
};

int main(){
    MyApp app;
    AudioDevice::printAll();

     //audioRate audioBlockSize audioOutputs audioInputs device

    //-1 for audioOutputs or audioInputs opens all available channels. see: AudioIO::channels()
    app.initAudio(SAMPLE_RATE, BLOCK_SIZE, 60, 0, -1);

    // Use this for 2809                   **CHANGE AUDIO INPUT TO 0 **
    //app.initAudio(SAMPLE_RATE, BLOCK_SIZE, 2, -1, AudioDevice("Aggregate Device").id());

    // Use this for sphere
    //    app.initAudio(44100, BLOCK_SIZE, -1, -1, AudioDevice("ECHO X5").id());

    app.start();
    return 0;
}
