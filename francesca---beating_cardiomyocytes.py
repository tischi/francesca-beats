from ij.io import OpenDialog
from ij.io import Opener
from ij.gui import GenericDialog
from ij.plugin import ZProjector, RGBStackMerge, SubstackMaker, Concatenator
from ij import IJ, ImagePlus, ImageStack
from ij.plugin import Duplicator
from ij.process import StackStatistics
from ij.plugin import ImageCalculator
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
import os, os.path, re, sys
from jarray import array
from ij.process import ImageConverter
import math
from math import sqrt
from ij.macro import MacroRunner

from loci.plugins import BF
from loci.common import Region
from loci.plugins.in import ImporterOptions

from automic.table import TableModel			# this class stores the data for the table
from automic.table import ManualControlFrame 	#this class visualises TableModel via GUI
from java.io import File
from ah.utils import ROIManipulator

# ParallelFFTJ
from edu.emory.mathcs.parallelfftj import FloatTransformer
from edu.emory.mathcs.parallelfftj import SpectrumType

# import my analysis function collection
import os, sys, inspect
this_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if this_folder not in sys.path:
  print this_folder
  sys.path.insert(0, this_folder)
import ct_analysis_functions as af
reload(af)

def mean(x):
  mean = sum(x) / len(x)
  return mean
  
def sd(x):
  mean = sum(x) / len(x)
  differences = [xx - mean for xx in x]
  sq_differences = [xx**2 for xx in differences]
  sd = sqrt(sum(sq_differences)/len(x))
  return sd

# Structure for batch analysis:
# 
# - main 
#  - parameters = get_analysis_parameters()
#  - folder = get_folder()
#  - table = init_results_table()
#  - data_info = get_data_info(folder) 
#  - batch_analyze(folder, parameters, table, dada) 
#    - iDataSet = 0
#    - for iDataSet in 'data_info' 
#       - imp, filename = return_imp(iDataSet, data_info) # might be different for different projects
#       - analyze(imp, folder, filename, parameters, table)
#         - do analysis on imp
#         - write results to table(iDataSet)
#         - write segmentation overlay images to filepath_output = f(folder, filename)
# - show interactive results table
# - save results table

# Improved structure:

# Structure for batch analysis:
#  
# - main 
#  - parameters = get_analysis_parameters()
#  - folder = get_folder()
#  - table = init_results_table()
#  - get_data_info(folder, table) 
#  - batch_analyze(parameters, table) 
#    - for row in table 
#       - analyze(row, table, parameters)
#         - imp = load_imp(table, row)
#         - write results to table(row)
#         - write segmentation overlay images (use from table(row))


def analyze(iDataSet, tbModel, p, output_folder):

  #
  # LOAD FILES
  #

  filepath = tbModel.getFileAPth(iDataSet, "RAW", "IMG")
  filename = tbModel.getFileName(iDataSet, "RAW", "IMG") 
  print("Loading: "+filepath)
  IJ.run("Bio-Formats Importer", "open=["+filepath+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
  imp = IJ.getImage()
  
  #
  # INIT
  #
  IJ.run("Options...", "iterations=1 count=1"); 

 
  #
  # SCALING
  #
  IJ.run(imp, "Scale...", "x="+str(p["scale"])+" y="+str(p["scale"])+" z=1.0 interpolation=Bilinear average process create"); 
  imp = IJ.getImage()
  # save output file
  output_file = filename+"--downscale_input.tif"
  IJ.saveAs(IJ.getImage(), "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "INPUT","IMG")

  #
  # CONVERSION
  #
  
  #IJ.run(imp, "8-bit", "");
 
  
  #
  # CROPPING
  #
  
  #imp.setRoi(392,386,750,762);
  #IJ.run(imp, "Crop", "");

  
  #
  # BACKGROUND SUBTRACTION
  #
  
  # IJ.run(imp, "Subtract...", "value=32768 stack");

  IJ.run(imp, "Z Project...", "projection=[Average Intensity]");
  imp_avg = IJ.getImage()
  ic = ImageCalculator();
  imp = ic.run("Subtract create 32-bit stack", imp, imp_avg);
 
  #
  # REGION SEGMENTATION
  #
  
  imp1 = Duplicator().run(imp, 1, imp.getImageStackSize()-1)
  imp2 = Duplicator().run(imp, 2, imp.getImageStackSize())
  imp_diff = ic.run("Subtract create 32-bit stack", imp1, imp2);
  #imp_diff.show()

  IJ.run(imp_diff, "Z Project...", "projection=[Standard Deviation]");
  imp_diff_sd = IJ.getImage()
 
  # save
  IJ.run(imp_diff_sd, "Gaussian Blur...", "sigma=5");
  output_file = filename+"--sd.tif"
  IJ.saveAs(imp_diff_sd, "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "SD","IMG")

  IJ.run(imp_diff_sd, "Enhance Contrast", "saturated=0.35");
  IJ.run(imp_diff_sd, "8-bit", "");
  IJ.run(imp_diff_sd, "Properties...", "unit=p pixel_width=1 pixel_height=1 voxel_depth=1");
  IJ.run(imp_diff_sd, "Auto Local Threshold", "method=Niblack radius=60 parameter_1=2 parameter_2=0 white");
  
  rm = ROIManipulator.getEmptyRm()
  IJ.run(imp_diff_sd, "Analyze Particles...", "add");


  # select N largest Rois
  diameter_roi = []
  for i in range(rm.getCount()):
    roi = rm.getRoi(i)
    diameter_roi.append([roi.getFeretsDiameter(), roi])
  diameter_roi = sorted(diameter_roi, reverse=True)
  #print diameter_roi

  rm.reset()
  for i in range(min(len(diameter_roi), p["n_rois"])):
    rm.addRoi(diameter_roi[i][1]) 
  
  # save 
  output_file = filename+"--rois"
  ROIManipulator.svRoisToFl(output_folder, output_file, rm.getRoisAsArray())  
  tbModel.setFileAPth(output_folder, output_file+".zip", iDataSet, "REGIONS","ROI")

   
  #
  # FFT in each region
  #

  IJ.run(imp, "Variance...", "radius=2 stack");
  output_file = filename+"--beats.tif"
  IJ.saveAs(imp, "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "BEATS","IMG")
  
  n = rm.getCount()
  for i_roi in range(n):
    imp_selection = Duplicator().run(imp)
    rm.select(imp_selection, i_roi)
    IJ.run(imp_selection, "Clear Outside", "stack");
    imp_selection.show()
    
    # FFT using Parallel FFTJ
    transformer = FloatTransformer(imp_selection.getStack())
    transformer.fft()
    imp_fft = transformer.toImagePlus(SpectrumType.FREQUENCY_SPECTRUM)
    imp_fft.show()

    # Analyze FFt
    IJ.run(imp_fft, "Gaussian Blur 3D...", "x=0 y=0 z=3");
    IJ.run(imp_fft, "Plot Z-axis Profile", "");
    output_file = filename+"--Region"+str(i_roi+1)+"--fft.tif"
    IJ.saveAs(IJ.getImage(), "TIFF", os.path.join(output_folder, output_file))
    tbModel.setFileAPth(output_folder, output_file, iDataSet, "FFT_R"+str(i_roi+1),"IMG")

    IJ.run(imp_fft, "Select All", "");
    rm.addRoi(imp_fft.getRoi())
    rm.select(rm.getCount())
    rt = ResultsTable()
    rt = rm.multiMeasure(imp_fft); #print(rt.getColumnHeadings);
    x = rt.getColumn(rt.getColumnIndex("Mean1"))
    #rm.runCommand("delete")
    
    peak_pos = []
    peak_height = []
    x_min = 10
    for i in range(x_min,len(x)/2):
      if (x[i]>x[i-1]) and (x[i]>x[i+1]):
        peak_pos.append(i)
        peak_height.append(float(x[i]))

    #print(peak_pos)
    #print(peak_height)
    peak_height, peak_pos = zip(*sorted(zip(peak_height, peak_pos), reverse=True))
    #print(peak_height)
    #print(peak_pos)
    #print(len(x))
    
    n_max = 3
    for i_max in range(min(len(peak_height),n_max)):
      tbModel.setNumVal(round(float(len(x))/float(peak_pos[i_max]),2), iDataSet, "F"+str(i_max+1)+"_R"+str(i_roi+1))
      tbModel.setNumVal(int(peak_height[i_max]), iDataSet, "A"+str(i_max+1)+"_R"+str(i_roi+1))
        

#
# ANALYZE INPUT FILES
#
def determine_input_files(foldername, tbModel):

  print("Determine input files in:",foldername)
  pattern = re.compile('(.*).czi') 
  #pattern = re.compile('(.*)--beats.tif') 
   
  i = 0
  for root, directories, filenames in os.walk(foldername):
	for filename in filenames:
	   print("Checking:", filename)
	   if filename == "Thumbs.db":
	     continue
	   match = re.search(pattern, filename)
	   if (match == None) or (match.group(1) == None):
	     continue
	   tbModel.addRow()
	   tbModel.setFileAPth(foldername, filename, i, "RAW","IMG")
	   print("Accepted:", filename)
	   
	   i += 1
    
  return(tbModel)
 
def get_parameters(p, num_data_sets):
  gd = GenericDialog("Correct 3D Drift Options")

  gd.addMessage("Found "+str(num_data_sets)+" data sets")
  gd.addStringField("analyse", "all");

  gd.addMessage("Image analysis parameters:")
  for k in p.keys():
    gd.addNumericField(k, p[k], 2);
  gd.showDialog()
  if gd.wasCanceled():
    return

  to_be_analyzed = gd.getNextString()
  for k in p.keys():
    p[k] = gd.getNextNumber()
    
  return to_be_analyzed, p

    
if __name__ == '__main__':

  print("")
  print("#")
  print("# Beating Analysis")
  print("#")
  print("")
  
  #
  # GET INPUT FOLDER
  #
  od = OpenDialog("Click on one of the image files in the folder to be analysed", None)
  input_folder = od.getDirectory()
  if input_folder is None:
    sys.exit("No folder selected!")
    

  #
  # MAKE OUTPUT FOLDER
  #
  output_folder = input_folder[:-1]+"--fiji"
  if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

  
  #
  # DETERMINE INPUT FILES
  #
  tbModel = TableModel(input_folder)
  tbModel.addFileColumns('RAW','IMG')
  tbModel = determine_input_files(input_folder, tbModel)

  #
  # GET PARAMETERS
  #
  p = dict()
  p["scale"] = 0.2
  p["n_rois"] = 2
 
  to_be_analyzed, p = get_parameters(p, tbModel.getRowCount())

  #
  # POPULATE AND SHOW INTERACTIVE TABLE
  #

  for i in range(p["n_rois"]):
    for j in range(3):
      tbModel.addValColumn("F"+str(j+1)+"_R"+str(i+1), "NUM")
      
  for i in range(p["n_rois"]):
    for j in range(3):
      tbModel.addValColumn("A"+str(j+1)+"_R"+str(i+1), "NUM")
 
  tbModel.addFileColumns('INPUT','IMG')
  tbModel.addFileColumns('BEATS','IMG')
  tbModel.addFileColumns('SD','IMG')
  tbModel.addFileColumns('REGIONS','ROI')
  
  for i in range(p["n_rois"]):
    tbModel.addFileColumns('FFT_R'+str(i+1),'IMG')
  
  frame=ManualControlFrame(tbModel)
  frame.setVisible(True)
  
  #
  # ANALYZE
  #

  if not to_be_analyzed=="all":
    af.close_all_image_windows()
    analyze(int(to_be_analyzed)-1, tbModel, p, output_folder)
  else:
    for i in range(tbModel.getRowCount()):
      af.close_all_image_windows()
      analyze(i, tbModel, p, output_folder)

  af.close_all_image_windows()