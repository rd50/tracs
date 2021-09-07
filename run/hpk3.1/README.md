## note
For running the simulation, just run: 

$ python config\_gene.py 

to start the related script.

But notice that the you will need a powerful server(~100GB memory and ~12 core) because large amount data will be calculated.

# To generate all the pictures we have for now, just follow these steps:

1. get the simulation results

$ cd ../../

$ source setup.sh

* set environment 

$ cd run/hpk3.1/

$ cd edge\_lgad\_1\_0\_changeV\_0000\_0um

* run HPK3.1 with CV doping simulation and wait for a relatively long time...

$ python config\_gene.py

$ cd ../edge\_lgad\_1\_0\_changeV\_4000\_2.2um

* run HPK3.1 with **modified** CV doping simulation and wait for a relatively long time...

$ python config\_gene.py

2. draw Vmax-z pictures at different bias voltage (z is the laser's vertical distance from sensor's top) 

$ python draw\_Vmax\_com.py

3. draw volt-time(waveform) pictures at 200V bias voltage

$ python draw\_wavf\_com.py

4. draw volt-time(waveform) pictures at 200V bias voltage

$ python draw\_chargecollection\_com.py

5. draw current contributions at 25um depth 200V bias voltage

$ root plt_current.c

6. draw electric field 

* You can set log Y-axis to make it more clear in interactive window.

$ root reverseEFX.c


