#pragma once

#include <juce_audio_utils/juce_audio_utils.h>
#include "PluginProcessor.h"

// Draws the big dark knob with the recessed pointer (design coordinates are 1280 x 1280).
class AngelLookAndFeel : public juce::LookAndFeel_V4
{
public:
    void drawRotarySlider (juce::Graphics&, int x, int y, int width, int height,
                           float sliderPosProportional, float rotaryStartAngle,
                           float rotaryEndAngle, juce::Slider&) override;
};

class AngelShiftEditor : public juce::AudioProcessorEditor,
                         private juce::Timer
{
public:
    explicit AngelShiftEditor (AngelShiftProcessor&);
    ~AngelShiftEditor() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    void timerCallback() override;

    static juce::Path makeWings();   // white outline, grows when the knob is turned right
    static juce::Path makeHorns();   // red, appears when the knob is turned left

    AngelShiftProcessor& processor;
    AngelLookAndFeel laf;
    juce::Slider knob;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> attachment;
    juce::Path wings, horns;
    double lastValue = 1.0e9;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (AngelShiftEditor)
};
