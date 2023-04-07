# LPC Vocoding
We generate LPC coefficients and gain estimates for windows of a given speech signal. We then vocode the speech signal by filtering an excitation with the LPC coefficients/gain for a frame. We introduce robustness by adding a 50% overlap between LPC frames so our output has twice the LPC coefficients used per frame.

To run, simply call main.m. It will compute the LPC coefficients for welcome16k.wav in 20ms frame intervals and generate two outputs. The first output is the the resynthesized signal vocoded with a 100Hz excitement. The other is the resynthesized signal vocoded with white noise.
