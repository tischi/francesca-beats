# francesca-beats

Quantification of beating cells using Fiji.

## Requirements

Fiji with java 1.8

## Installation

- dowload and unzip: https://github.com/tischi/francesca-beats/archive/master.zip
- copy the .jar file to your Fiji plugins folder 
- copy the folder named to your Fiji plugins folder
- copy the two .py files together in one arbiratry folder

## Usage

- start Fiji
- drag&drop the francesca*.py file onto Fiji and click Run at the bottom of the script editor

## Output

The script finds the N (input parameter __n rois__) largest region centers where beating occurs. 
It measures the beating frequencies in these region centers.

The output values are as such:

- F1_R1_NUM = frequency 1 in region 1; to facilitate comparison with raw data the units are __period length in frames__; to compute actual frequencies in Hz you have to compute 1/(F1_R1_NUM*dt), where dt is the time difference between frames.
- A1_R1_NUM = amplitude of frequency 1 

Frequencies are sorted according to their amplitude => F1 has the highest amplitude

Regions are sorted according to their diameter => R1 has the largest diameter





