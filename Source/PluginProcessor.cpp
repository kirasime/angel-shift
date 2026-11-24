#include "PluginProcessor.h"
#include "PluginEditor.h"

AngelShiftProcessor::AngelShiftProcessor()
    : AudioProcessor (BusesProperties()
                          .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                          .withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      apvts (*this, nullptr, "STATE", createLayout())
{
    shiftParam = apvts.getRawParameterValue (kShiftID);
    setLatencySamples (0);   // zero-latency by design
}

juce::AudioProcessorValueTreeState::ParameterLayout AngelShiftProcessor::createLayout()
{
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

    params.push_back (std::make_unique<juce::AudioParameterFloat> (
        juce::ParameterID { kShiftID, 1 },
        "Formant Shift",
        juce::NormalisableRange<float> (-12.0f, 12.0f, 0.01f),
        0.0f,
        juce::AudioParameterFloatAttributes().withLabel ("st")));

    return { params.begin(), params.end() };
}

void AngelShiftProcessor::prepareToPlay (double sampleRate, int)
{
    shifter.setShiftSemitones (shiftParam->load());
    shifter.prepare (sampleRate);
}

bool AngelShiftProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    return layouts.getMainInputChannelSet()  == juce::AudioChannelSet::stereo()
        && layouts.getMainOutputChannelSet() == juce::AudioChannelSet::stereo();
}

void AngelShiftProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;

    const int numSamples  = buffer.getNumSamples();
    const int numChannels = buffer.getNumChannels();

    for (int ch = 2; ch < numChannels; ++ch)
        buffer.clear (ch, 0, numSamples);

    if (numChannels < 1)
        return;

    shifter.setShiftSemitones (shiftParam->load());
    shifter.process (buffer.getWritePointer (0),
                     numChannels > 1 ? buffer.getWritePointer (1) : nullptr,
                     numSamples);
}

juce::AudioProcessorEditor* AngelShiftProcessor::createEditor()
{
    return new AngelShiftEditor (*this);
}

void AngelShiftProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    if (auto xml = apvts.copyState().createXml())
        copyXmlToBinary (*xml, destData);
}

void AngelShiftProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    if (auto xml = getXmlFromBinary (data, sizeInBytes))
        if (xml->hasTagName (apvts.state.getType()))
            apvts.replaceState (juce::ValueTree::fromXml (*xml));
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new AngelShiftProcessor();
}
