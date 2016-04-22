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
import os
import os.path
import re
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

# import my analysis function collection
import os, sys, inspect
this_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if this_folder not in sys.path:
  print this_folder
  sys.path.insert(0, this_folder)
import ct_analysis_functions as af
reload(af)

all_results = []

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
  
  IJ.run(imp, "8-bit", "");
 
  
  #
  # CROPPING
  #
  
  #imp.setRoi(392,386,750,762);
  #IJ.run(imp, "Crop", "");

  
  #
  # BACKGROUND SUBTRACTION
  #
  
  # IJ.run(imp, "Subtract...", "value=32768 stack");
  
  #
  # ENHANCEMENT
  #
  
  IJ.run(imp, "Stack Difference", "gap=1"); imp = IJ.getImage()
  IJ.run(imp, "32-bit", "");
  IJ.run(imp, "Variance...", "radius=5 stack");

  # make output file
  output_file = filename+"--beats.tif"
  IJ.saveAs(IJ.getImage(), "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "BEATS","IMG")
    
  #
  # SEGMENTATION
  # 

  #IJ.setMinAndMax(imp, 0, 4095); IJ.run(imp, "16-bit", "");
  #IJ.run(imp, "3D Spot Segmentation", "seeds_threshold="+str(p["threshold"])+" local_background=0 radius_0=2 radius_1=20 radius_2=22 weigth=0.5 radius_max=10 sd_value=1 local_threshold=[Local Mean] seg_spot=Block volume_min=1 volume_max=1000000 seeds=Automatic spots="+imp.getTitle()+"radius_for_seeds=2 output=[Label Image]");
  # get number of objects in 3D
  #imp_label = IJ.getImage()
  #num_objects = StackStatistics(imp_label).max
  #imp_bw = af.threshold(imp_label, 1)

  #
  # MEASURE
  #

  IJ.run(imp, "Plot Z-axis Profile", "");
  output_file = filename+"--plot.tif"
  IJ.saveAs(IJ.getImage(), "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "PLOT","IMG")

  IJ.run(imp, "Select All", "");
  rm = RoiManager(True)
  rm.addRoi(imp.getRoi());
  rt = rm.multiMeasure(imp); #print(rt.getColumnHeadings);
  x = rt.getColumn(rt.getColumnIndex("Mean1"))

  peak_pos = []
  peak_height = []
  for i in range(1,len(x)-1):
    if (x[i]>x[i-1]) and (x[i]>x[i+1]):
      peak_pos.append(i)
      peak_height.append(float(x[i]))
  peak_dt = []
  for i in range(1,len(peak_pos)):
    peak_dt.append(float(peak_pos[i]-peak_pos[i-1]))

  def mean(x):
    mean = sum(x) / len(x)
    return mean
  
  def sd(x):
    mean = sum(x) / len(x)
    differences = [xx - mean for xx in x]
    sq_differences = [xx**2 for xx in differences]
    sd = sqrt(sum(sq_differences)/len(x))
    return sd

  tbModel.setNumVal(round(mean(peak_dt),3), iDataSet, "PeakDiff_mean")
  tbModel.setNumVal(round(sd(peak_dt),3), iDataSet, "PeakDiff_sd")
  tbModel.setNumVal(round(mean(peak_height),3), iDataSet, "PeakHeight_mean")
  tbModel.setNumVal(round(sd(peak_height),3), iDataSet, "PeakHeight_sd")
  
  #
  # EVALUATE MEASUREMENTS
  #

  #print rt.getColumnHeadings()
  
  '''
  issues = ""
  tbModel.setNumVal(num_objects, iDataSet, "Segmented_Particles")
  
  if num_objects == 0:
    issues = issues + " no_particles "  
  elif num_objects > 1:
    issues = issues + " multiple_particles " 
  else:    
    # not round
    rc  = rt.getColumn(rt.getColumnIndex("Round"))
    tbModel.setNumVal(min(rc), iDataSet, "Minimal_Roundness")
    if min(rc) < p["minimal_roundness"]:
      issues = issues + " elongated "      
    # shifting
    xy  = zip(rt.getColumn(rt.getColumnIndex("X")), rt.getColumn(rt.getColumnIndex("Y")))
    shifts = [math.sqrt((v[0]-xy[0][0])**2 + (v[1]-xy[0][1])**2) for v in xy]
    tbModel.setNumVal(max(shifts), iDataSet, "Max_Shift")
    if max(shifts) > p["maximal_shift"]:
      issues = issues + " shifting "
 
  tbModel.setValue(issues, iDataSet, "Issues")
  tbModel.setBoolVal(len(issues)==0, iDataSet, "OK")
  '''
  

#
# ANALYZE INPUT FILES
#
def determine_input_files(foldername, tbModel):
  
  pattern = re.compile('(.*).czi')  
  i = 0
  for root, directories, filenames in os.walk(foldername):
	for filename in filenames:
	   print(filename)
	   if filename == "Thumbs.db":
	     continue
	   match = re.search(pattern, filename)
	   if (match == None) or (match.group(1) == None):
	     continue
	   tbModel.addRow()
	   tbModel.setFileAPth(foldername, filename, i, "RAW","IMG")
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
  
  #
  # GET INPUT FOLDER
  #
  od = OpenDialog("Click on one of the image files in the folder to be analysed", None	)
  input_folder = od.getDirectory()
  #input_folder = "C:/Users/almf/Desktop/kaia/data/"
  if input_folder is None:
    fff

  #
  # MAKE OUTPUT FOLDER
  #
  output_folder = input_folder[:-1]+"--fiji"
  if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

 
  #
  # INIT OUTPUT TABLE
  #
  tbModel = TableModel(input_folder)
  tbModel.addFileColumns('RAW','IMG')
  tbModel.addFileColumns('INPUT','IMG')
  tbModel.addFileColumns('BEATS','IMG')
  tbModel.addFileColumns('PLOT','IMG')
  tbModel.addValColumn("PeakDiff_mean", "NUM")
  tbModel.addValColumn("PeakDiff_sd", "NUM")
  tbModel.addValColumn("PeakHeight_mean", "NUM")
  tbModel.addValColumn("PeakHeight_sd", "NUM") 
  tbModel.addColumn("Issues")
  tbModel.addValColumn("OK", "BOOL")
  frame=ManualControlFrame(tbModel)
  frame.setVisible(True)

  #
  # DETERMINE INPUT FILES
  #
  tbModel = determine_input_files(input_folder, tbModel)
  
  #
  # GET PARAMETERS
  #
  p = dict()
  p["scale"] = 0.2
    
  to_be_analyzed, p = get_parameters(p, tbModel.getRowCount())
  
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