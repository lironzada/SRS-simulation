README SRS Imaging SIMULATION
Written by Bart Fokker, Vrije Universiteit Amsterdam, 21 February 2019
Edited by Liron Zada, December 2019

Purpose: Simulating an SRS setup to check settings before starting the measurements and to investigating the results
of a setup type before building it.

Description: Simulation which can simulate up to 6 SRS channels. The simulation can be 2 dimensional with cylinders
representing beads or 1 dimensional with planes representing single particles. Here, we simulate an SRS setup that
extracts the SRS from the noise through a Lock-in Amplifier. The low pass filter that is simulated is a digital RC
filter, since our setup uses the HF2LI of Zurich Instruments.

Files needed: 	ImagingSimulationInterface.mlapp
		SimulationLIA.m
		Stokes modulation schematic.PNG
		Pump and Stokes modulation schematic.PNG
		Double Stokes modulation schematic.PNG

For interface, open file: ImagingSimulationInterface.mlapp

Control parameters:

    Setup Tab

      - Modulation method: Single Carrier encoding, Double Carrier encoding and Probe and Carrier encoding. The former
	is the common method, the Carrier beam is modulated and the Probe is measured. The second method has two 
	modulations on the Carrier beam, 1 with a high frequency, the other with several low frequencies so that the
	different channels are encoded in the sidebands. The last is similar to the second case, only with the Probe
	beam producing the sidebands for the different channels.

      - Sample choice: 2D bead, 1D particle.

      - Phase encoding: single or double phase encoding. For double phase encoding two channels are encoded at the
	same frequency, but 90 degrees out of phase. These channels do not interact, so with the same total 
	frequency range, more distance between the channels becomes possible, decreasing the crosstalk.

      - Data collection: digital data collection (1 sample) or analog data collection (averaging). Our setup has 2
	data collection methods: the first method sends the LIA output directly to the computer sending only 1 LIA 
	output sample per pixel. The second method creates an analog signal, which is send to an ADC, which averages
	all the LIA samples for each pixel.

    Input parameters Tab

      - Laser powers: average Probe and Carrier laser power on the sample in mW

      - Pixel size: Size of pixels

      - Number of pixels: number of x and y pixels

      - Pixel Dwell TIme: pixel dwell time

      - Filter order: filter order of the cascaded digital RC filter

      - Time Constant: Time constant of the digital RC filter

    Setup parameters Tab

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
	The focussed Gaussian beam model is used in this simulation. w0x and w0y are the beam waists in the x and
	y direction respectively. The beamwidth is defined as: Intensity = exp(-2 X^2 / w^2), with w as the width,
	which is equal to the waist at the focus. Multiplying the waists of both the Carrier and Probe beams,
	gives the effective SRS beam waist, which is shown in the interface.

    Channels Tab

      - Number of channels: Number of channels limited from 1 to 6.

      - Modulation frequency: modulation frequency of the first channel

      - Channel distance: distance between channels in Hz

      - Sample properties: Contains a list to which samples can be added an deleted. The properties of each sample
        can be set with the 4 fields at the right side of the list. These fields are:

         - Channel: Specifies to which channel the sample belongs.

         - Radius or object size: Radius of the bead (2D sample) or the object size (1D sample).

         - X or Y position: Gives the bead centre (2D sample) or the offset of the object (1D sample).

    Advanced Settings Tab

      - Advanced settings switch: This switch controlls if the advanced settings are used, or the ones from the
        other tabs.

      - Number of channels: Number of channels limited from 1 to 6.

      - Table containing the Csignal matrix: This table enables the user to simulate an SRS setup more closely by
        having materials being active in more than 1 demodulation channel. The rows represent the wavenumbers of
        the input channels (laser channels). The columns represent the sample types (materials), so 1 column is a
        selection of the Raman spectrum of that material.

      - Carrier/Probe Power: Here the power of each laser channel can be set independently. So if the modulation
        method is the Single Carrier modulation or de Double Carrier modulation, it controlls the carrier powers.
        If it is Probe and Carrier modulation, the fields control the Probe powers. The power of the non-split
        laser is still controlled by the usual field.

      - Mod frequency: Modulation frequency of each (laser) channel. If a double modulation method is applied,
        these fields become the low frequency components, which form channels through the sidebands with the fast
        modulation frequency.

      - Fast modulation frequency: This field is only visible when a double modulation method is selected. It
        controls the fast modulation frequency.

      - Phase: This field controls the phase of the modulation and demodulation of the channels. These can obly be
        set to 0 (in-phase) and pi/2 (quadrature). These fields are only selectable when double phase encoding is
        selected.

Input support:

      - Simulation image preview: For visualizing the sample before starting the simulation. This way, it is easier
	to get the desired test sample. Time can be saved this way.

      - Size: Size of image in um, cannot be changed, is only meant as support for the input settings.

Output: Simulation image

Other info:

	The simulation can be started by pushing on the 'Start Simulation' button.

	Simulation progress: The first part of the simulation simply shows the message 'Calculating...'. When the
	demodulation part of the simulation is reached, the message shows the percentage of how much is
	demodulated already.

	Simulation time: To get an indication of how long the simulation can take, we mention here our simulation
	times. For 30x30 pixels, a PDT of 100 microseconds, 6 channels, a filter order of 8 and a pixel size of 1
	micrometer, our simulation takes 70 seconds. For the same settings but with 100x100 pixels, it takes us
	770 seconds.