# -------------------------------------------------------------
# Script Name: Ratiometric Image Conversion Script
# Version: 16
# Description:
#     This script ratiometrically converts confocal images 
#     acquired in 2 channels.
# -------------------------------------------------------------

#@ String(value="<html style=\"width: 400px;text-align: center;\"><p style=\"margin:0px;padding:0px;font-size:12px;\"><b>Ratiometric image conversion in Fiji</b></p><p>This script ratiometrically converts confocal images acquired in 2 channels.<br/>Hover over (?) for more information.<br/><br/></p></html>", visibility=MESSAGE, required=false) desc
#@ File(label="Folder with your images (?)", style="directory", description="Input folder") src_dir
#@ File(label="Folder to save your images (?)", style="directory", description="Output folder") out_dir
#@ String(label="Extension for the images to look for (?)", value="tif", description="Extension of your images to select in the input folder") extension
#@ String(value="<html style=\"width: 240px;\"><div style=\"height: 1px; background: #c0c0c0;\"/></html>", visibility=MESSAGE, required=false) line1
#@ Integer(label="Channel to use for segmentation (?)", value=1, description="For best results select the channel using 405 nm laser", min=1, max=2) seg_chnl
#@ String (label="Channel ratio (?)", choices={"ch2/ch1", "ch1/ch2"}, description="Decide on which channel should be divided by which. For HPTS image ratio divide the 458 nm by the 405 nm channel") ratio_method
#@ Boolean(label="Do channel alignment (?)", value=False, description="<html>Enable for the recursive alignment of a stack of images<br/><br/><b>Requires to run:</b><ol><li><b>TurboReg</b> (https://bigwww.epfl.ch/thevenaz/turboreg/)</li><li><b>StackReg</b> (https://bigwww.epfl.ch/thevenaz/stackreg/)</li></ol></html>") do_stackreg
#@ String(value="<html style=\"width: 240px;\"><div style=\"height: 1px; background: #c0c0c0;\"/></html>", visibility=MESSAGE, required=false) line2
#@ String (label="How to deal with thresholding (?)", choices={"Fully automatic", "Manual once and apply to all", "Fully manual"}, description="<html>Select a mode of action to set thresholds for the images<br/><br/><b>Requires to run:</b><ol><li><b>ImageScience</b> (https://imagej.net/libs/imagescience)</li></ol></html>") thresh_method
#@ String (label="Select LUT (?)", choices={"Green Fire Blue", "Fire", "Grays", "Ice", "Spectrum", "Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "Red/Green", "Cyan Hot", "HiLo", "ICA", "ICA2", "ICA3", "Magenta Hot", "Orange Hot", "Rainbow RGB", "Red Hot", "Thermal", "Yellow Hot", "blue orange icb", "cool", "gem", "glow", "mpl-inferno", "mpl-magma", "mpl-plasma", "mpl-viridis", "phase", "physics", "royal", "sepia", "smart", "thal", "thallium", "unionjack"}, description="Select the lookup table for final coloring") lut_method
#@ String (label="Image Processing Range Max (?)", choices={"Default and apply to all", "Manual once and apply to all", "Fully manual"}, description="<html>Select a mode of action to set the max value used to set the <br/>highest color in the LUT<br/><ul><li>Default values: min=0 and max=3</li></ol></html>") image_processing
#@ String(value="<html style=\"width: 400px;text-align: center;\">Please cite Barbez et al. 2017<br/>and Rößling et al. 2025</html>", visibility=MESSAGE, required=false) footer

# ─── IMPORTS ────────────────────────────────────────────────────────────────────

import os
import sys
import time
import fnmatch

from java.lang import Double

from ij import IJ, Menus, Prefs, WindowManager
from ij.gui import WaitForUserDialog, HTMLDialog, GenericDialog
from ij.plugin import Duplicator, ImageCalculator, LutLoader
from ij.plugin.frame import ThresholdAdjuster

from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader
from loci.formats import MetadataTools


# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────

def getFileList(directory, extensions):
    """Get a list of files with the extension

    Parameters
    ----------
    directory : str
        Path of the files to look at
    extensions : [str]
        Extensions to look for

    Returns
    -------
    list
        List of files with the extension in the folder
    """
    
    files = []
    extensions = [ext.lower() for ext in extensions]
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            # Check if the file matches any of the extensions (case-insensitive)
            if any(fnmatch.fnmatch(filename.lower(), "*" + ext) for ext in extensions):
                files.append(os.path.join(dirpath, filename))
    return files


def progress_bar(progress, total, line_number, prefix=""):
    """Progress bar for the IJ log window

    Parameters
    ----------
    progress : int
        Current step of the loop
    total : int
        Total number of steps for the loop
    line_number : int
        Number of the line to be updated
    prefix : str, optional
        Text to use before the progress bar, by default ''
    """

    size = 30
    x = int(size * progress / total)

    IJ.log(
        "\\Update%i:%s\t[%s%s] %i/%i\r"
        % (line_number, prefix, "#" * x, "." * (size - x), progress, total)
    )


def get_series_info_from_ome_metadata(path_to_file):
    """Get the number of series from a file

    Parameters
    ----------
    path_to_file : str
        Path to the file

    Returns
    -------
    int
        Number of series for the file
    """
    reader = ImageReader()
    reader.setFlattenedResolutions(False)
    omeMeta = MetadataTools.createOMEXMLMetadata()
    reader.setMetadataStore(omeMeta)
    reader.setId(path_to_file)
    series_count = reader.getSeriesCount()

    series_index = []
    for i in range(series_count):
        if i == 0:
            resolution_count = 0
            series_index.append(resolution_count)
        else:
            reader.setSeries(i - 1)
            resolution_count += reader.getResolutionCount()
            series_index.append(resolution_count)

    reader.close()

    return series_count, series_index


def open_single_series_with_BF(path_to_file, series_number):
    """Open a single serie for a file using Bio-Formats

    Parameters
    ----------
    path_to_file : str
        Path to the file
    series_number : int
        Number of the serie to open

    Returns
    -------
    ImagePlus
        ImagePlus of the serie
    """
    options = ImporterOptions()
    options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE)
    options.setSeriesOn(series_number, True)
    
    options.setId(path_to_file)
    imps = BF.openImagePlus(options)  # is an array of imp with one entry
    return imps[0]


def timed_log(message):
    """Print a log message with a timestamp added

    Parameters
    ----------
    message : str
        Message to print
    """
    IJ.log(time.strftime("%H:%M:%S", time.localtime()) + ": " + message)
    
    
def command_exists(c):
    """Check if a command exists. If not, the according plugin may need to be installed

    Parameters
    ----------
    c : str
        Name of the command

    Returns
    -------
    Boolean
        True if command exists
    """
    return c in Menus.getCommands()
    

def close_window(name):
    """Closes a window that may be open

    Parameters
    ----------
    name : str
        Name of the window
    """
    window_thres = WindowManager.getWindow("Threshold")
    if window_thres is not None:
        window_thres.close()



# ─── Variables ────────────────────────────────────────────────────────────────

# Predefine the amount of Gaussian Blur
blur_radius = 1.5
# Predefine the sigma for the Laplace Filter (higher value for bigger structures)
laplace_radius = 1
# Predefine the Radius for Background subtraction
back_subtract = 20
# Predefine the threshold method
threshold_method = "Default"

# Minimum and Maximum Value for output Image
min_display_default = 0
max_display_default = Double(3)
min_display = min_display_default
max_display = max_display_default
manually_set_display = False

# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

IJ.log("\\Clear")
timed_log("Script starting")

src_dir = src_dir.getAbsolutePath()
out_dir = out_dir.getAbsolutePath()

filename_filter = [extension]
files = getFileList(src_dir, filename_filter)

thresh_value = 0


# Check if user wants auto alignment but lacks the StackReg/TurboReg plugin
# ! StackReg & TurboReg command names contain a space at the end
if do_stackreg and (not command_exists("StackReg ") or not command_exists("TurboReg ")):
	HTMLDialog("Stackreg Plugin Missing", "<html><body style=\"width:350px;\">To enable channel alignment of a stack of images, the StackReg plugin needs to be installed. The StackReg plugin requires the installation of a second plugin called TurboReg. Otherwise, deactivate the channel alignment option and run the script again.<ol><li><b>TurboReg</b> (<a href=\"https://bigwww.epfl.ch/thevenaz/turboreg/\">https://bigwww.epfl.ch/thevenaz/turboreg/</a>)</li><li><b>StackReg</b> (<a href=\"https://bigwww.epfl.ch/thevenaz/stackreg/\">https://bigwww.epfl.ch/thevenaz/stackreg/</a>)</li></ol></body></html>")

# Check if user has ImageScience Plugin installed
elif not command_exists("FeatureJ Laplacian"):
	HTMLDialog("ImageScience Plugin Missing", "<html><body style=\"width:350px;\">For the thresholding to work, the ImageScience plugin needs to be installed. Install the plugin and run the script again.<ol><li><b>ImageScience</b> (<a href=\"https://imagej.net/libs/imagescience\">https://imagej.net/libs/imagescience</a>)</ol></body></html>")

# If the list of files is not empty
elif files:

    # Check if selected LUT exists. Otherwise use default "Green Fire Blue"
    if not LutLoader.getLut(lut_method):
    	# Check if Green Fire Blue exists 
        if LutLoader.getLut("Green Fire Blue"):
            timed_log("LUT '" + lut_method + "' does not exist. Using default: 'Green Fire Blue'")
            lut_method = "Green Fire Blue"
        else:
        	# Green Fire Blue does not exist. Use a built-in LUT
            timed_log("LUT '" + lut_method + "' does not exist. Using default: 'Fire'")
            lut_method = "Fire"


    for file_id, file in enumerate(sorted(files)):
        # Get info for the files
        folder = os.path.dirname(file)
        basename = os.path.basename(file)
        
        progress_bar(file_id + 1, len(files), 1, "Processing: " + str(file_id))

        # Import the file with BioFormats
        series_count, series_index = get_series_info_from_ome_metadata(file)
        for series in range(series_count):
            progress_bar(series + 1, series_count, 2, "Opening series : ")

            imp = open_single_series_with_BF(file, series_index[series])
            
            # Get the number of z-slices
            n_slices = imp.getNSlices()
            for z in range(1, n_slices + 1):
                progress_bar(z, n_slices, 3, "Opening z-slice in series : ")

                # Set the Z slice
                imp.setZ(z)

	            # If the imp only has one channel, we avoid operations that require at least two channels
                is_single_channel = imp.getNChannels() < 2
                if is_single_channel:
	                seg_chnl = 1

                # Align channels if needed
                if do_stackreg and not is_single_channel:
	                out_ratio = os.path.join(out_dir, imp.getTitle() + "_z" + str(z) + "_ratio_chnlaligned.tif")
	                IJ.run(imp, "StackReg ", "transformation=[Affine]")
                else:
	                out_ratio = os.path.join(out_dir, imp.getTitle() + "_z" + str(z) + "_ratio.tif")

                IJ.run(imp, "32-bit", "")
                IJ.run(imp, "Gaussian Blur...", "sigma=" + str(blur_radius) + " stack")

                imp_c1 = Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
                imp_c1.setTitle("C1")

                imp_c2 = Duplicator().run(imp, 2, 2, 1, 1, 1, 1)
                imp_c2.setTitle("C2")

                imp_segment = Duplicator().run(imp, seg_chnl, seg_chnl, 1, 1, 1, 1)
                imp_segment.setTitle("Backsubtract")
                imp_segment.show()
                IJ.selectWindow(imp_segment.getID()) 
                IJ.run(imp_segment, "Subtract Background...", "rolling=" + str(back_subtract))
                IJ.run("FeatureJ Laplacian", "compute smoothing=" + str(Double(laplace_radius)))
                imp_laplace = IJ.getImage()

                # Hide Backsubtract image
                imp_segment.changes = False
                imp_segment.hide()


                # Ask for thresholds
                if thresh_method == "Fully automatic":
	            	# Hide Laplace image in automatic mode
	            	imp_laplace.hide()

	            	# Apply auto threshold to the Laplace image
	                IJ.setAutoThreshold(imp_laplace, threshold_method)

	                # Get threshold values from the Laplace image
	                thresh_value = imp_laplace.getProcessor().getMaxThreshold()
	            	min_thresh_value = imp_laplace.getProcessor().getMinThreshold()
                elif thresh_value == 0:
	            	# Pre-set thresholds.
	            	IJ.setAutoThreshold(imp_laplace, threshold_method)

	                # Open Threshold window and bring it to the front
	                IJ.run("Threshold...")
	                IJ.selectWindow("Threshold")
	                
	                # Pre-define settings in the Threshold window
	            	ThresholdAdjuster.setMode("Red")
	            	ThresholdAdjuster.setMethod(threshold_method)
	                
	                # Wait for user to set manual threshold...
	                WaitForUserDialog("Set manual threshold", "Select the \"Threshold\" window that opened in the background\nand set the according threshold. Press OK to continue the script\n*after* setting the threshold.").show()
	                
	                # Get the manually set threshold values from the Laplace image
	                thresh_value = imp_laplace.getProcessor().getMaxThreshold()
	            	min_thresh_value = imp_laplace.getProcessor().getMinThreshold()
	                
	                # Hide the Laplace image
	                imp_laplace.hide()
	                
	
				# Set Threshold in case of "Manual once and apply to all"
                IJ.setThreshold(imp_laplace, min_thresh_value, thresh_value)
                
                # Reset threshold for the next image
                if thresh_method == "Fully manual":
	                thresh_value = 0
                else:
	                # Close threshold window
	                close_window("Threshold")
	            		
	            
	            # Ask for image_processing values  
                ask_for_image_processing = image_processing == "Manual once and apply to all" or image_processing == "Fully manual"
                if image_processing == "Manual once and apply to all" and manually_set_display:
	            	# Image processing values have already been asked
	            	ask_for_image_processing = False
                if ask_for_image_processing:
	            	manually_set_display = True
	            	
	            	# Reset back to default values
	            	min_display = min_display_default
	            	max_display = max_display_default
	            	
	            	# Wait for user to set image_processing values
	            	gui = GenericDialog("Image Processing Range")
	
	            	# Add some gui elements (Ok and Cancel button are present by default)
	            	# Elements are stacked on top of each others by default (unless specified)
	            	gui.setInsets(0,0,0)
	            	if image_processing == "Fully manual":
	            		gui.addMessage("Set the max value for '" + imp.getTitle() + "'\nthat are used to set highest color in the chosen LUT:")
	            	else:
	            		gui.addMessage("Set the max value to set the highest\ncolor in the chosen LUT that are used for all images:")
	            	
	            	gui.setInsets(5,10,0)
	            	gui.addNumericField("max", max_display_default) 
	            	gui.hideCancelButton()
	            	gui.showDialog()
	            	
	            	if gui.wasOKed():
	                    # Get the selected max value
	                    max_display = gui.getNextNumber()


	            # Convert the image to black and white. The mask will have an inverting 
                # LUT (white is 0 and black is 255) unless Black Background is checked
                Prefs.blackBackground = False
                IJ.run(imp_laplace, "Convert to Mask", "")
                if imp_laplace.isInvertedLut():
	                IJ.run(imp_laplace, "Invert LUT", "")


                IJ.run(imp_laplace, "32-bit", "")
                IJ.setAutoThreshold(imp_laplace, "Default dark")
                IJ.run(imp_laplace, "NaN Background", "")
                IJ.run(imp_laplace, "Divide...", "value=255")


				# Create Ratiometric Image
				# Divide the signal intensity of the 458-nm channel by the 405-nm channel
                if ratio_method == "ch2/ch1":
                    imp_result = ImageCalculator().run("Divide create 32-bit", imp_c2, imp_c1)
                elif ratio_method == "ch1/ch2":
	                imp_result = ImageCalculator().run("Divide create 32-bit", imp_c1, imp_c2)
                imp_result = ImageCalculator().run("Multiply create 32-bit", imp_result, imp_laplace)

				# Set Display and save
                IJ.run(imp_result, lut_method, "")
                IJ.run(imp_result, "Select None", "")
                IJ.setMinAndMax(imp_result, min_display, max_display)

                # Set calibration bar
                IJ.run(
	                imp_result,
	                "Calibration Bar...",
	                "location=[Upper Right] fill=White label=Black number=5 decimal=3 font=12 zoom=1 overlay",
	            )
                IJ.saveAs(imp_result, "Tiff", out_ratio)


	         	# Close all images
                imp.changes = False
                imp.close()
                imp_result.changes = False
                imp_result.close()
                imp_c1.changes = False
                imp_c1.close()
                imp_c2.changes = False
                imp_c2.close()
                imp_laplace.changes = False
                imp_laplace.close()
                imp_segment.changes = False
                imp_segment.close()


# Close threshold window
close_window("Threshold")

timed_log("Script finished !")
print("Script finished !")
