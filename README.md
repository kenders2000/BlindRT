Blind Estimation of Reveberation Time from speech and music
===========================================
Paul kendrick, University of Salford, Acoustics research centre, for more information please contact p.kendrick@salford.ac.uk

This program can compute the reverberation time from passivley recieved signals, such as speech and music. The advantage of this methodology is that it can be performed in-situ, with the room occupied. This code is free to use, under the condition that a reference is made to:

*Paul Kendrick, Nicola Shiers, Robert Conetta, Trevor J. Cox, Bridget M. Shield, Charlie Mydlarz*, Blind estimation of reverberation time in classrooms and hospital wards, Applied Acoustics, Volume 73, Issue 8, August 2012, Pages 770-780, ISSN 0003-682X, http://dx.doi.org/10.1016/j.apacoust.2012.02.010.


Description of algorithm
--------------
The algorithm is written in matlab, it contains a number of subfunctions, the following list describes the operation of the algorithm, for more details please consult the references at the end of this readme.
  - Choose wav file, this needs to be a recording of speech or music,
  for speech at least an hour, for music perhaps longer
 -  read in section of wav file, and filter into octave band
 -     Maximum likehood estimation of decays; 
         -- Polyfit algorithm, this uses 0.5 s windows on the envelope of
          the signal, fitting a 1st order poly nomial to the log
          envelope, the gradient of each fit is extracted, windows are
          moved along, using a very high overlap (0.002s), and regions
          which are continuously decaying identified from the gradient
          (gradient equivlent for RT==100s used as maximum gradient) -
          these decay phases are then fine tuned (start maximum, end min
          of LP filtered envelope)
         -- Maximum Likelihood fit to every decay phase, a model of sound
         decay in a room is where the envelope, alpha*a.^I+(1-alpha)*b.^I,
         modulates gaussian white noise.  This is a
         model of non-diffuse sound decay, utilising two exponential
         decays, with decay rates determined by a and b, added together
         using a convex sum where the balence between the two is
         controlled by the parameter alpha, all parmaters are constrained
         between 0 & 1. The fit proceedure first computes the likelihood
         over a coarse grid of values for a and b (alpha optimsed for
         each grid point), then performs a fine search using the corase
         result as a starting point, the function fmincon is used.  This
         is contrained minimisation.
         --   Once all ML parameters are computed, the dynamic range of
         each decayphase is computed from the ML decay curve.
  -     Postprocess results, for every octave band, a framework for
  estimationg the RT is as follows, find the length of signal required to
  ensure that at least 40 decay phases with at least 25 dB of dynamic
  range are present, (this is just a rule of thumb), I found for speech
  3-5 (segLen is this parameter - in secs) mins works.  In each 3-5 min 
 section, compute the RT from  the 
 decay phases with at least 25dB dynamic range, select the minimum RT as
 the estimate for that section.  Compute the minimum RT for multiple
 sections.  Average over all sections provides the blind RT estimate, the
 95% CL from the standard error indicates the confidence in that result.
 
 References
 ------
-*Paul Kendrick, Nicola Shiers, Robert Conetta, Trevor J. Cox, Bridget M. Shield, Charlie Mydlarz*, Blind estimation of reverberation time in classrooms and hospital wards, Applied Acoustics, Volume 73, Issue 8, August 2012, Pages 770-780, ISSN 0003-682X, http://dx.doi.org/10.1016/j.apacoust.2012.02.010.
-*P Kendrick, TJ Cox, FF Li, Y Zhang, JA Chambers*, Monaural room acoustic parameters from music and speech - The Journal of the Acoustical Society of America, 2008
-*P Kendrick, FF Li, TJ Cox, Y Zhang, JA Chambers*, Blind estimation of reverberation parameters for non-diffuse rooms, Acta Acustica united with Acustica, 2007
