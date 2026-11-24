#include "PluginEditor.h"

namespace
{
    // Design space: 1280 x 1280 (the size of the reference pictures).
    constexpr float kDesign   = 1280.0f;
    constexpr float kKnobCx   = 640.0f;
    constexpr float kKnobCy   = 660.0f;
    constexpr float kKnobR    = 282.0f;
    constexpr float kStroke   = 22.0f;

    const juce::Colour kKnobOuter  { 0xff202020 };
    const juce::Colour kKnobInner  { 0xff2a2a2a };
    const juce::Colour kTitle      { 0xffcdeeff };
    const juce::Colour kHorn       { 0xffff0000 };
}

//==============================================================================
void AngelLookAndFeel::drawRotarySlider (juce::Graphics& g, int x, int y, int width, int height,
                                         float pos, float startAngle, float endAngle, juce::Slider&)
{
    const auto bounds = juce::Rectangle<float> ((float) x, (float) y, (float) width, (float) height);
    const auto c = bounds.getCentre();
    const float R = juce::jmin (bounds.getWidth(), bounds.getHeight()) * 0.5f;
    const float k = R / kKnobR;

    g.setColour (kKnobOuter);
    g.fillEllipse (c.x - R, c.y - R, 2.0f * R, 2.0f * R);

    const float Ri = R - 24.0f * k;
    g.setColour (kKnobInner);
    g.fillEllipse (c.x - Ri, c.y - Ri, 2.0f * Ri, 2.0f * Ri);

    // recessed pointer: hollow rectangle with a dark outline, pointing up at 0 st
    juce::Path pointer;
    pointer.addRectangle (-62.0f * k, -255.0f * k, 124.0f * k, 215.0f * k);
    const float angle = startAngle + pos * (endAngle - startAngle);
    pointer.applyTransform (juce::AffineTransform::rotation (angle).translated (c.x, c.y));

    g.setColour (kKnobInner);
    g.fillPath (pointer);
    g.setColour (kKnobOuter.darker (0.15f));
    g.strokePath (pointer, juce::PathStrokeType (22.0f * k, juce::PathStrokeType::mitered));
}

//==============================================================================
AngelShiftEditor::AngelShiftEditor (AngelShiftProcessor& p)
    : juce::AudioProcessorEditor (&p), processor (p)
{
    wings = makeWings();
    horns = makeHorns();

    knob.setLookAndFeel (&laf);
    knob.setSliderStyle (juce::Slider::RotaryVerticalDrag);
    knob.setTextBoxStyle (juce::Slider::NoTextBox, false, 0, 0);
    knob.setRotaryParameters (juce::MathConstants<float>::pi * 1.25f,
                              juce::MathConstants<float>::pi * 2.75f, true);
    knob.setDoubleClickReturnValue (true, 0.0);
    knob.setMouseDragSensitivity (300);
    knob.setPopupDisplayEnabled (true, true, this);
    knob.setTextValueSuffix (" st");
    addAndMakeVisible (knob);

    attachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment> (
        processor.apvts, AngelShiftProcessor::kShiftID, knob);

    setResizable (false, false);
    setSize (600, 600);
    startTimerHz (30);
}

AngelShiftEditor::~AngelShiftEditor()
{
    stopTimer();
    knob.setLookAndFeel (nullptr);
}

void AngelShiftEditor::timerCallback()
{
    // host automation moves the slider without notifying us: repaint the wings / horns
    if (knob.getValue() != lastValue)
    {
        lastValue = knob.getValue();
        repaint();
    }
}

void AngelShiftEditor::resized()
{
    const float s = (float) getWidth() / kDesign;
    knob.setBounds (juce::Rectangle<float> ((kKnobCx - kKnobR) * s, (kKnobCy - kKnobR) * s,
                                            2.0f * kKnobR * s, 2.0f * kKnobR * s).toNearestInt());
}

//==============================================================================
juce::Path AngelShiftEditor::makeWings()
{
    // left wing, traced by eye from the reference picture (design coordinates)
    juce::Path left;
    left.startNewSubPath (395.0f, 505.0f);
    left.cubicTo (300.0f, 455.0f, 235.0f, 330.0f, 222.0f, 190.0f);
    left.cubicTo (200.0f, 330.0f, 190.0f, 450.0f, 255.0f, 578.0f);
    left.lineTo  (150.0f, 568.0f);
    left.lineTo  (128.0f, 495.0f);
    left.cubicTo (170.0f, 640.0f, 205.0f, 720.0f, 205.0f, 790.0f);
    left.lineTo  (65.0f, 745.0f);
    left.cubicTo (80.0f, 900.0f, 180.0f, 1010.0f, 215.0f, 1050.0f);
    left.lineTo  (280.0f, 1135.0f);
    left.cubicTo (330.0f, 1150.0f, 390.0f, 1110.0f, 425.0f, 1075.0f);

    juce::Path both (left);
    left.applyTransform (juce::AffineTransform::scale (-1.0f, 1.0f).translated (kDesign, 0.0f));   // mirror
    both.addPath (left);
    return both;
}

juce::Path AngelShiftEditor::makeHorns()
{
    juce::Path left;
    left.startNewSubPath (185.0f, 300.0f);                                  // outer curve
    left.cubicTo (120.0f, 420.0f, 150.0f, 540.0f, 305.0f, 605.0f);
    left.startNewSubPath (185.0f, 300.0f);                                  // inner curve
    left.cubicTo (190.0f, 420.0f, 250.0f, 490.0f, 355.0f, 518.0f);

    juce::Path both (left);
    left.applyTransform (juce::AffineTransform::scale (-1.0f, 1.0f).translated (kDesign, 0.0f));
    both.addPath (left);
    return both;
}

void AngelShiftEditor::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colours::black);

    const float s = (float) getWidth() / kDesign;

    // title
    g.setColour (kTitle);
    g.setFont (juce::Font (juce::FontOptions ("Segoe Script", 128.0f * s, juce::Font::plain)));
    g.drawText ("Angel Shift", juce::Rectangle<float> (0.0f, 120.0f * s, (float) getWidth(), 220.0f * s),
                juce::Justification::centred, false);

    // wings (right turn = formants up) / horns (left turn = formants down)
    const float v = (float) knob.getValue();
    const float amount = juce::jlimit (0.0f, 1.0f, std::abs (v) / 12.0f);
    if (amount < 0.01f)
        return;

    const float shown = std::pow (amount, 0.6f);         // appears quickly, saturates towards the end
    const float grow  = 0.88f + 0.12f * shown;           // slight growth around the knob
    const auto tf = juce::AffineTransform::scale (grow, grow, kKnobCx, kKnobCy).scaled (s);

    juce::Path p (v > 0.0f ? wings : horns);
    p.applyTransform (tf);
    g.setColour ((v > 0.0f ? juce::Colours::white : kHorn).withAlpha (shown));
    g.strokePath (p, juce::PathStrokeType (kStroke * s * grow, juce::PathStrokeType::curved,
                                           juce::PathStrokeType::rounded));
}
