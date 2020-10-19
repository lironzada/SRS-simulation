README OPTIMAL SETTINGS
Written by Bart Fokker, Vrije Universiteit Amsterdam, 21 February 2019

Purpose: Obtain the optimal settings for SRS setup

Description: 	In this simulation, the optimal settings are calculated by finding the settings for which the CNR
		is maximal while still retaining as much responsivity as possible. The smaller the filter 
		bandwidth is, the less noise will be transmitted, increasing the CNR. However, the smaller 
		bandwidth results in a slower response, which may cause the loss of details/resolution. Therefore,
		an equilibirum needs to be found for responsivity and noise suppression. Here, we defined the 
		optimal setting as the maximal CNR. The CNR is the SNR times the maximum-minimum aka 
		constrast devided by sqrt(2) to account the two noise terms (1 for max and the other for min).
		This simulation assumes a digital RC filter, since our LIA (the HF2LI) has such a filter.

Files needed: 	OptimalSettingsInterface.mlapp
		SystematicSettings.m
		LIAOptimalSettings.m

For interface, open file: OptimalSettingsInterface.mlapp

Control parameters:

    Parameter Tab 1
	
      - 2 optimization choices:

		CNR or Imaging Time: Define the CNR or the Time
			When choosing Time: Define Imaging Time or Pixel Dwell Time (both in seconds)

		Tightly packed or single particles: optimize for CNR
		

      - Flyback: Signifies how long the flyback is in % of the Pixel Dwell Time

      - xpixels: number of x pixels

      - ypixels: number of y pixels

      - pixelsize: the size of 1 pixel in micrometer

      - objectsize: the length of the object for which the TC is optimized. Also called the aimed resolution.
   
      -- xFOV and yFOV: Field of view, which is an 'output' included for clarity and cannot be directly changed.
	

    Parameter Tab 2

      - filter order: selects the filter order of the cascaded digital RC filter

      - Find optimal filter order: to find the optimal filter order, this box can be checked. Then, a new field 
	appears	after the filter order field. By entering a number (B) in this field that is higher than the filter 
	order (A), the simulation is run with several filter orders (from A to B). From these, the optimal filter 
	order is selected by choosing the one with the largest Nett SNR.

      - total power: The total power sets the total laser power that hits the sample in mW. There is a choice to set the
	Carrier and Probe powers individually or to set the total power on the sample. The total power will then 
	divide the power over the 2 laser beams with a ratio of 1:2 (Probe:Carrier) since this results in the largest
	SNR.

    Setup parameter Tab

      - Coefficients
	Csignal, Cshot and Cthermal are coefficients we measured for our setup for Polystyrene beads.
	They are defined as:

		signal = Csignal x ProbePower x CarrierPower
		shot noise = Cshot x sqrt(ProbePower x NEPBW)
		thermal noise = Cthermal x sqrt(NEPBW)

		with NEPBW as the Noise Equivalent Power Bandwidth, which is the bandwidth of a block filter that
		transmits just as much noise as the actual filter. Multiplying the square root of the NEPBW with 
		the Spectral Density gives the noise that the filter transmits.
		NEPBW = Gamma(order - 1/2) / (4 x sqrt(pi) x TC x Gamma(order))

      - Beam parameters
	In this simulation, the focussed Gaussian beam model is used. w01 and w02 are the beam waist for beam 1 and
	2 respectively. The beamwidth is defined as: Intensity = exp(-2 X^2 / w^2), with w as the width, which is
	equal to the beam waist in the focus.


Outputs:
	
	Imaging Time: 	time to collect 1 frame
	SNR: 		signal to noise ratio
	Nett SNR:	signal to noise ratio multiplied with the normalized peak of peak-valley value
	Filter order: 	Optimal filter order
	TC: 		Time constant of the filter
	PDT: 		Pixel dwell time
	Carrier power: 	Carrier power on the sample
	Probe power: 	Probe power on the sample

	Also a figure is plotted showing of the simulation, to illustrate the results.