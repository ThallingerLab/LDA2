/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class TooltipTexts
{
  
  public final static String TABS_MAIN_QUANTITATION = "<html>Quantitation of a single experiment</html>";
  public final static String TABS_MAIN_BATCH = "<html>Batch Quantitation of a experiments</html>";
  public final static String TABS_MAIN_STATISTICS = "<html>Analysis of quantified results</html>";
  public final static String TABS_MAIN_DISPLAY = "<html>Manual verification of results with 3D viewer</html>";
  public final static String TABS_MAIN_SETTINGS = "<html>Default LDA settings; e.g. type of mass spectrometer</html>";
  public final static String TABS_MAIN_LICENSE = "<html>Information about the current licensing status</html>";
  public final static String TABS_MAIN_HELP = "<html>Help component</html>";
  public final static String TABS_MAIN_ABOUT = "<html>About</html>";
  
  public final static String TABS_RESULTS_SELECTION = "<html>Select here the files you want to analyze</html>";
  public final static String TABS_RESULTS_GROUP = "<html>Heat maps and bar charts for ";
  public final static String TABS_RESULTS_OVERVIEW = "<html>Lipid classes overview</html>";
  public final static String TABS_RESULTS_OVERVIEW_GROUPS= "<html>Lipid classes overview for the selected groups</html>";
  public final static String TABS_RESULTS_HEATMAP = "<html>Heat map of ";
  public final static String TABS_RESULTS_BARCHART = "<html>Bar chart of ";
  public final static String TABS_RESULTS_HEATMAP_GROUP = "<html>Heat map of groups of ";
  public final static String TABS_RESULTS_BARCHART_GROUP = "<html>Bar chart of groups of ";
  
  public final static String ACCEPT_GENERAL = "<html>Accept settings</html>";
  public final static String CANCEL_GENERAL = "<html>Cancel operation</html>";
  public final static String GENERAL_ACCEPT_SINGLE_START = "<html>Accept ";
  public final static String GENERAL_ACCEPT_SINGLE_END = "?</html>";
  public final static String ALL_GENERAL  = "<html>Select all</html>";
  public final static String NONE_GENERAL  = "<html>Select none</html>";
  public final static String INVERT_GENERAL  = "<html>Invert selection</html>";;
  
  public final static String QUANTITATION_SINGLE_RAW_FILE = "<html>Enter the raw file of your experiment!<br>" +
  		  "The following file types are accepted:<br>" +
  		  "RAW-file:   XCalibur RAW files (XCalibur must be installed)<br>" +
  		  "RAW-dir:    MassLynx RAW directory (MassLynx must be installed)<br>" +
  		  "mzXML:      mzXML format<br>" +
  		  "chrom:      LDA quant file format</html>"; 
  
  public final static String QUANTITATION_SINGLE_MASS_LIST ="<html>Enter the mass list file containing the analytes of interest!<br>" +
    "The file type must be Microsoft Excel (xls or xlsx).<br>" +
    "A detailed description how to prepare such a file can be found in the<br>" +
    "user manual Appendix A or in the Help component.<br>" +
    "Examples are provided in the \"examples\" directory of the installation</html>";
  
  public final static String QUANTITATION_RET_BEFORE = "<html>Retention time tolerance in minutes.<br>"+
    "Specification of a time tolerance to allow quantitation<br/>"+
    "before the expected retention time entered in the Excel file.<br/>" +
    "If no time is entered in the Excel, this field has no effect!</html>";
  
  public final static String QUANTITATION_RET_AFTER = "<html>Retention time tolerance in minutes.<br>"+
  "Specification of a time tolerance to allow quantitation<br/>"+
  "after the expected retention time entered in the Excel file.<br/>" +
  "If no time is entered in the Excel, this field has no effect!</html>";
  
  public final static String QUANTITATION_CUTOFF = "<html>Peaks smaller than \"highest peak area\" times this <br>"+
  "value in per mille are automatically discared.</html>";
  
  public final static String QUANTIFICATION_RET_SHIFT = "<html>Retention time shift in minutes relative to the one<br>"+
  "entered in the mass list Excel file</html>";
  
  public final static String QUANTITATION_ISOTOPES = "<html>How many isotopic peaks should be quantified additional to the base peak.<br/>" +
  		"The first input field defines the amount of isotpes.<br/>" +
  		"The second one the amount of isotopes that must match the<br/>" +
  		"isotopic distribution (otherwise hit is discarded).</html>";
  
  public final static String QUANTITATION_PROCESSORS = "<html>The number of processors to be used for quantitation.<br/>" +
      "It is recommended to use n-1 processors.<br/>" +
      "In this way the remaining system can run properly.</html>";
  
  public final static String QUANTITATION_ION_MODE = "<html>The used ion mode for this search.<br/>" +
      "Alex123 target lists contain both ion modes, thus,<br/>" +
      "the one for the search has to be selected.</html>";

  
  public final static String QUANTITATION_RET_UNKNOWN = "<html>Specifies wether or not a molecule should be quantified where no <br/>" +
  "retention time has been entered in the mass list Excel file.</html>";
  
  public final static String QUANTITATION_START = "<html>Starts the quantitation.<br/>" +
  "Before the quantitation can start a RAW and an Excel file has to be selected!</html>";
  
  public final static String QUANTITATION_STATUS_TEXT = "Current quantitation status";
  public final static String QUANTITATION_PROGRESS = "Quantitation progress";
  
  public final static String QUANTITATION_BATCH_RAW_FILE = "<html>Enter a folder containing the raw files of your experiment!<br>" +
  "The following file types are accepted as input:<br>" +
  "RAW-file:   XCalibur RAW files (XCalibur must be installed)<<br>" +
  "RAW-dir:    MassLynx RAW directory (MassLynx must be installed)<br>" +
  "mzXML:      mzXML format<br>" +
  "chrom:      LDA quant file format</html>";
  
  public final static String QUANTITATION_BATCH_MASS_LIST = "<html>Enter a folder containing the mass list file with the analytes of interest!<br>" +
  "The file type must be Microsoft Excel (xls or xlsx).<br>" +
  "A detailed description how to prepare such a file can be found in the<br>" +
  "user manual Appendix A or in the Help component.<br>" +
  "Examples are provided in the \"examples\" directory of the installation</html>";
  
  public final static String QUANTITATION_BATCH_START = "<html>Starts the quantitation.<br/>" +
  "Before the quantitation can start a RAW and an Excel files have to be selected!</html>";
  
  public final static String STATISTICS_ADD_RESULT_FILES = "<html>Add files for the results analysis.<br>" +
  "Single or multiple result files from the quantitation output can be selected.<br>" +
  "The file type must be Microsoft Excel (xls or xlsx).<br>" +
  "The files are added in the table below</html>";
  
  public final static String STATISTICS_ADD_DIR = "<html>Add files for the results analysis.<br/>" +
  "A directory containing the result files from the quantitation has to be selected.<br/>" +
  "The program reads the files in the foldert automatically (not from subfolders).<br/>" +
  "The file type must be Microsoft Excel (xls or xlsx).<br/>" +
  "The files are added in the table below</html>";
  
  public final static String STATISTICS_REMOVE_SELECTION =  "<html>Removes selected files from the table below.<br/>" +
  "The items have to be selected before in the table below.</html>";
  
  public final static String STATISTICS_REMOVE_ALL_SELECTION =  "<html>Removes all files from the table below.<br/>" +
  "It cleans the list.</html>";
  
  public final static String STATISTICS_ADD_TO_GROUP =  "<html>Combines the selected files to a group.<br/>" +
  "This is necessary to get statistics about experiments belonging together.</html>";
  
  public final static String STATISTICS_GROUP_SELECTION_NAME = "<html>The name of the group.</html>";
  public final static String STATISTICS_GROUP_SELECTION_RENAME = "<html>Rename the group.</html>";
  public final static String STATISTICS_GROUP_SELECTION_REMOVE = "<html>Rename the group.</html>";
  public final static String STATISTICS_GROUP_SELECTION_DELETE = "<html>Deletes the enitre group.</html>";
  
  public final static String STATISTICS_SEPARATE_RT = "<html>The algorithm can find several peaks for one analyte<br/>" +
  "With this option the user can specify if these are displayed as separate hits or not.<br/>" +
  "This option has no effects for data was generated with LDA versions 1.6.0 or earlier.<br/>" +
  "The value specifies the RT-shift-tolerance of the peak grouping accross experiments.</html>";

  
  public final static String STATISTICS_CORRECT_ORDER_FILE = "<html>The analytes are normally ordered according to the order in the result Excel files.<br/>" +
  "However, due to not-quantified analytes this order could deviate from the original order in the quantitation file.<br/>" +
  "In order to guarantee the same order, the original quantitation file can be specified here.<br/>" +
  "This value is not mandatory!</html>";
  public final static String STATISTICS_IS_PREFIX = "<html>The prefix for the internal standard (usage of standards is not mandatory).<br/>" +
  "The standards in the result Excel file should carry a certain prefix in order to determine<br/>" +
  "which result is a standard and which not.</html>";
  
  public final static String STATISTICS_ES_PREFIX = "<html>The prefix for the external standard (usage of standards is not mandatory).<br/>" +
  "The standards in the result Excel file should carry a certain prefix in order to determine<br/>" +
  "which result is a standard and which not.</html>";
  
  public final static String STATISTICS_ADD_ABS_SETTINGS = "<html>Settings for absolute standard concentrations, lipid or protein content (not mandatory).<br/>" +
  "These settings allow standardization of the data on different values (standards, protein content, etc.).</html>";
  public final static String STATISTICS_ADD_CUTOFF_SETTINGS = "<html>Settings for class specific cutoff thresholds (not mandatory).<br/>" +
  "For each lipid class, small intensities can be filtered out.<br/>" +
  "The threshold must be exceeded in one sample only, then, the analyte remains in all samples</html>";
  public final static String STATISTICS_REMOVE_ABS_SETTINGS = "<html>Removes the entered absolute settings.</html>";
  public final static String STATISTICS_REMOVE_CUTOFF_SETTINGS = "<html>Removes the entered cutoff settings.</html>";
  
  public final static String STATISTICS_ACCEPT_SELECTION = "<html>Accept selection and calculate statistics.</html>";
  
  public final static String STATISTICS_ABS_LOAD = "<html>Loads previously stored settings (wqs.xml).</html>";
  public final static String STATISTICS_ABS_SAVE = "<html>Saves entered absolute settings for later reuse.</html>";

  public final static String STATISTICS_CUTOFF_LOAD = "<html>Loads previously stored settings (cut.xml).</html>";
  public final static String STATISTICS_CUTOFF_SAVE = "<html>Saves entered absolute settings for later reuse.</html>";
  
  public final static String STATISTICS_ABS_TAB_CLASS = "<html>Standard settings for ";
  public final static String STATISTICS_ABS_TAB_EXP= "<html>Volume and concentration settings for experiment ";
  
  public final static String STATISTICS_ABS_STANDARDS_EXPS_SAME = "<html>Are the concentrations and volumes of the standards the same in all experiments?</html>";
  public final static String STATISTICS_ABS_STANDARDS_STDS_SAME = "<html>Do all of the standards have the same volume and concentration?</html>";

  public final static String STATISTICS_ABS_STANDARDS_TAB_EXP= "<html>Standard settings for experiment ";
  
  public final static String STATISTICS_ABS_STANDARDS_DILLUTION = "<html>Dillution from the time point where the external standard has been added to the internal standard</br>" +
  		"This option is additionally useful for the comparison of differntly dilluted experiments.</html>";
  
  public final static String STATISTICS_ABS_APPLY_TO_ALL = "<html>Apply this value to all input fields of this kind.</html>";
  public final static String STATISTICS_ABS_APPLY_TO_GROUP = "<html>Apply this value to all group members this experiment belongs to.</html>";
  
  public final static String STATISTICS_ABS_EXT_VOLUME = "<html>Volume of external standard</html>";
  public final static String STATISTICS_ABS_EXT_CONC = "<html>Concentration of external standard</html>";
  
  public final static String STATISTICS_ABS_INT_VOLUME = "<html>Volume of internal standard</html>";
  public final static String STATISTICS_ABS_INT_CONC = "<html>Concentration of internal standard</html>";
  
  public final static String STATISTICS_ABS_SAMPLE_VOLUME = "<html>The extracted volume before further processing (time point of adding of external standard).</html>";
  public final static String STATISTICS_ABS_SAMPLE_WEIGHT = "<html>The extracted sample weight (losses of analytes before addition of external standard are neglected).</html>";
  public final static String STATISTICS_ABS_END_VOLUME = "<html>The extracted volume after processing (time point of adding of internal standard).</html>";
  public final static String STATISTICS_ABS_PROTEIN_CONC = "<html>The protein content in the sample.</html>";
  public final static String STATISTICS_ABS_LIPID_CONC = "<html>The neutral lipid content in the sample.</html>";
  
  public final static String EXPORT_GENERAL = "<html>The available export options.</html>";
  public final static String EXPORT_PNG = "<html>Export image to PNG format (portable image format).</html>";
  public final static String EXPORT_SVG = "<html>Export image to SVG format (scalable vector graphics).</html>";
  public final static String EXPORT_EXCEL = "<html>Export data in Excel format.</html>";
  public final static String EXPORT_MZTAB = "<html>Export data in mzTab format.</html>";
  public final static String EXPORT_RDB = "<html>Export data in relational database format.</html>";
  public final static String EXPORT_TXT = "<html>Export data tab-delimited plain text format.</html>";
  public final static String EXPORT_CHROMS = "<html>Export the chromatograms in a picture.</html>";
  public final static String EXPORT_N_RTs = "<html>Export mass list for &omega identification.</html>";
  public final static String EXPORT_MAF = "<html>Exports data in MAF format.</html>";
  
  public final static String HEATMAP_SHOW_INT = "<html>Display the internal standards in the heat map or not.</html>";
  public final static String HEATMAP_SHOW_EXT = "<html>Display the external standards in the heat map or not.</html>";
  public final static String HEATMAP_ISOTOPES = "<html>Upper limit for the additional isotopic peaks used for the quantitative value.</html>";
  public final static String HEATMAP_DOUBLE_PEAKS = "<html>Highlight double peaks with a yellow rectangle.</html>";
  
  public final static String HEATMAP_SETTINGS = "<html>Settings of how to calculate the quantitative value.</html>";
  public final static String HEATMAP_SELECTED = "<html>Selection of molecules that should be displayed in the heat map.</html>";
  public final static String HEATMAP_COMBINED = "<html>Selection of molecules to create a combined bar chart.</html>";
  public final static String HEATMAP_EXPORT_OPTIONS = "<html>Options for the data export (which values should be exported, what should be in the raw and what in the column).</html>";

  public final static String HEATMAP_VALUE_TYPE = "<html>The quantitative type that is used for the geration of the heat map.<br/>" +
  "A detailed description about the types can be found in the user<br>" +
  "manual chapter 5.2  or in the Help component of the software.</html>";
  
  public final static String HEATMAP_IS_CORRECTION = "<html>Should the values be corrected with the internal standard?</html>";
  public final static String HEATMAP_ES_CORRECTION = "<html>Should the values be corrected with the external standard?</html>";
  public final static String HEATMAP_STANDARD_CORRECTION_TYPE = "<html>The possible correction types are:<br/>" +
  		"most reliable standard: correction on LDA method<br/>" +
  		"median: correction on the median of all available standars<br/>" +
  		"correction on a single standard</html>";
  public final static String HEATMAP_CONSIDER_DILUTION =  "<html>Should the dilution be considered?</html>";
  public final static String HEATMAP_DIVISOR_UNIT =  "<html>The display units for the divisors; e.g. m for /mL of a concentration</html>";
  public final static String HEATMAP_USE_AU =  "<html>Should arbitrary units be used or absolute entities?</html>";
  
  public final static String HEATMAP_EXPORT_DEVIATION = "<html>Should the deviation values be exported (standard deviation or standar error)?</html>";
  public final static String HEATMAP_EXPORT_STANDARD_DEVIATION = "<html>Use the standard deviation?</html>";
  public final static String HEATMAP_EXPORT_STANDARD_DEVIATION_WHICH = "<html>Single standard deviation, double standard deviation etc?</html>";
  public final static String HEATMAP_EXPORT_STANDARD_ERROR = "<html>Use the standard error?</html>";
  public final static String HEATMAP_EXPORT_COLUMN_ANALYTE = "<html>Table setup: analytes in columns and experiments/groups in rows</html>";
  public final static String HEATMAP_EXPORT_COLUMN_EXPERIMENT = "<html>Table setup: experiments/groups in columns and analytes in rows</html>";
  public final static String HEATMAP_EXPORT_RT = "<html>Should the retention time be exported?</html>";
  public final static String HEATMAP_EXPORT_RT_SD = "<html>Should the standard deviation of the retention time be exported?</html>";
  public final static String HEATMAP_EXPORT_SPECIES = "<html>The exported shall be on the species level</html>";
  public final static String HEATMAP_EXPORT_CHAIN = "<html>The exported shall be on the chain level</html>";
  public final static String HEATMAP_EXPORT_POSITION = "<html>The exported shall be on the position level.<br/>"
      + "A position will be exported if a majority for a position is available<br/>"
      + "and no contradicting evidence.</html>";
  
  public final static String BARCHART_VALUE_TYPE = "<html>The quantitative type that is used for the geration of the bar chart.<br/>" +
  "A detailed description about the types can be found in the user<br>" +
  "manual chapter 5.3  or in the Help component of the software.</html>";
  public final static String BARCHART_MAGNITUDE_CHOOSER = "<html>Selection of corresponding unit magnifier.</html>";
  public final static String BARCHART_SINGLE_SIDED = "<html>Paints the bar chart into one (positive) direction</html>";
  public final static String BARCHART_DOUBLE_SIDED = "<html>Paints the bar chart relative to the reference value (pos and neg direction).</html>";
  public final static String BARCHART_LOGARITHMIC = "<html>Paints the values in logarithmic scale.</html>";
  public final static String BARCHART_LINEAR = "<html>Paints the values in linear scale.</html>";
  public final static String BARCHART_LOG10 = "<html>Paints the values in decadic logarithm.</html>";
  public final static String BARCHART_LOG2 = "<html>Paints the values in binary scale.</html>";
  public final static String BARCHART_STANDARD_DEVIATION = "<html>Paints standard deviation as error bars.</html>";
  public final static String BARCHART_STANDARD_DEVIATION_WHICH = "<html>Multiplicative factor for the error bars of the standard deviation.</html>";
  public final static String BARCHART_STANDARD_ERROR = "<html>Paints standard error as error bars.</html>";
  
  public final static String OVERVIEW_CONSIDER_STANDARD = "<html>Compare the analyte classes based on standards, or on the area without standardization</html>";
  public final static String OVERVIEW_CONSIDER_DILUTION =  "<html>Should the dilution be considered for the comparison of the classes?</html>";
  public final static String OVERVIEW_CHOOSE_SELECTED = "<html>Selection of experiments/groups that should be displayed in the heat map.</html>";
  
  public final static String DISPLAY_OPEN_CHROM = "<html>Enter the raw file of your experiment! Just chrom file is allowed</html>";
  public final static String DISPLAY_OPEN_RESULT = "<html>Enter the result file of the quantitation! Just Excel file is allowed</html>";
  public final static String DISPLAY_START = "<html>Starts the display! Chrom and result file have to be selected before!</html>";
  public final static String DISPLAY_SELECT_CLASS = "<html>Selecte the analyte class you want to inspect.</html>";
  public final static String DISPLAY_SELECTION_TABLE = "<html>Click left mouse button to select an analyte and display it.<br/>" +
  		"Click right mouse button to add analyte before or after, or delete it.</html>";
  public final static String DISPLAY_MZ_PLUS = "<html>Display range for the 3D viewer in Dalton.<br/>" +
  		"It defines which m/z range after the m/z of the analyte should be displayed</html>";
  public final static String DISPLAY_MZ_MINUS = "<html>Display range for the 3D viewer in Dalton.<br/>" +
  "It defines which m/z range before the m/z of the analyte should be displayed</html>";
  public final static String DISPLAY_UPDATE = "<html>Update the settings in the 3D viewer.</html>";
  public final static String DISPLAY_SHOW_2D = "<html>Should the 2D viewer be displayed or not?<br/>" +
  		"The 2D viewer allows to change the quantitation manually. The 3D viewer is just for displaying purposes</html>";
  public final static String DISPLAY_SHOW_MSN = "<html>Should the names be displayed using MSn evidence?</html>";
  public final static String DISPLAY_LOCK_MZ = "<html>This setting locks the 3D viewer to a certain m/z range</html>";
  public final static String DISPLAY_RT_START = "<html>Display range for the 3D viewer in minutes.<br/>" +
  "This defines the start retention time.</html>";
  public final static String DISPLAY_RT_STOP = "<html>Display range for the 3D viewer in minutes.<br/>" +
  "This defines the stop retention time.</html>";


  public final static String DISPLAY_ADD_ANALYTE_NAME = "<html>The name of the analyte (plus double bonds if appropriate).</html>";
  public final static String DISPLAY_ADD_ANALYTE_FORMULA = "<html>The chemical formula of the analyte (without modification)<br/>" +
  		"The formula has to be entered by the chemical symbol<br/>" +
  		"followed by the cardinality, next chemical symbol and so on.</html>"; 
  public final static String DISPLAY_ADD_ANALYTE_MOD_NAME = "<html>The name of the modification.</html>";
  public final static String DISPLAY_ADD_ANALYTE_MOD_FORMULA = "<html>The chemical formula of the modification.<br/>" +
  		"For the reduction of elements place a \"-\" in front of the group. E.g., \"-NH4\" <br/>" +
  		"means the loss of one N and of 4 H, which would be the same as \"-N -H4\",<br/>" +
  		" wheras \"-N H4\" would mean the loss of one N and the addition of 4H.</html>";
  public final static String DISPLAY_ADD_ANALYTE_MZ = "<html>The m/z width of the primary chromatogram.<br/>" +
  "In order to be comparable with the other quantitations, it is not advised to change this value</html>";
  public final static String DISPLAY_ADD_ANALYTE_CHARGE = "<html>The charge of the analyte.</html>";
  public final static String DISPLAY_ADD_ANALYTE_RT = "<html>The RT identifier of the analyte.</html>";
  public final static String DISPLAY_ADD_ANALYTE_MASS = "<html>The m/z of the analyte.</html>";
  public final static String DISPLAY_ADD_ANALYTE_ACCEPT = "<html>Accept the settings and add the analyte.</html>";
  public final static String DISPLAY_ADD_ANALYTE_OH = "<html>The number of hydroxylation sites.</html>";
  
  public final static String DIALOG_ADD_PEAK_START_TIME = "<html>The lower time border in minutes.</html>";
  public final static String DIALOG_ADD_PEAK_STOP_TIME = "<html>The upper time border in minutes.</html>";
  public final static String DIALOG_ADD_PEAK_START_MZ= "<html>The lower m/z border in Da.</html>";
  public final static String DIALOG_ADD_PEAK_STOP_MZ= "<html>The upper m/z border in Da.</html>";
  public final static String DIALOG_ADD_PEAK_3D= "<html>Should the peak be defined in 3D (m/z borders settable)?</html>";
  
  public final static String SETTINGS_MS_MACHINE = "<html>The machine type you want to use for your quantitation</html>";
  public final static String SETTINGS_BUTTON_APPLY = "<html>Apply the selected settings to the current session</html>";
  public final static String SETTINGS_BUTTON_SAVE = "<html>Save the selected settings; after restart, these settings will be loaded</html>";
  public final static String SETTINGS_MSN_FRAGMENTATION = "<html>The MSn fragmentation settings you want to use for your quantitation</html>";

  
  public final static String HELP_HELP_OPEN= "<html>This link opens the help component.</html>";
  public final static String HELP_USER_MANUAL= "<html>These links return the user manual.</html>";
  public final static String HELP_EXAMPLES= "<html>These links return the documents that describe the LDA with the aid of example data.</html>";
  public final static String HELP_EXAMPLE_DOWNLOAD= "<html>Link to the download of example data.<br/>" +
  		"ATTENTION: Not all of the data is necessary to get an idea how the LDA works!<br/>" +
  		"Take a look at the examples document (previous link) to avoid unnecessary downloads.</html>";
  
  public final static String TABS_CHROMEXPORT_EXPS = "<html>Select samples to export</html>";
  public final static String TABS_CHROMEXPORT_ANALS = "<html>Select analytes to export</html>";
  public final static String TABS_CHROMEXPORT_MODS = "<html>Select modifications to export</html>";

  public final static String EXPORT_STATUS_TEXT = "<html>Current export status</html>";
  public final static String EXPORT_PROGRESS = "<html>Quantitation progress</html>";
  public final static String EXPORT_STOP = "<html>Interrupts the currently running export<html/>";
  
  public final static String RDI_GENERAL_SINGLECHAIN = "<html>Is an identification on one found chain sufficient?<html/>";
  public final static String RDI_GENERAL_CHAINCUTOFF = "<html>Cutoff for removing minor chains<br/>The value is relative to the sum of intensities of the highest identified chain combination<html/>";
  public final static String RDI_GENERAL_RTMAXDEV = "<html>Maximally allowed time deviation from the calcuated RT model<br/>If nothing is entered, the algorithm automatically allows four times the mean deviation as tolerance<html/>";
  public final static String RDI_GENERAL_MSIDORDER = "<html>The sequence for the application of the MS1 and MSn algorithms (influences search time):<br/>"
      + "&nbsp;&nbsp;MS<sup>1</sup> first: the data is searched first with the MS<sup>1</sup> algorithm first - default setting<br/>"
      + "&nbsp;&nbsp;MS<sup>n</sup> only: MSn spectra are searched for candidates - MS<sup>1</sup> identification is performed only where MS<sup>n</sup> hits are detected<br/>"
      + "&nbsp;&nbsp;MS<sup>n</sup> first: MSn spectra are searched for candidates - MS<sup>1</sup> identidication is performed where MS<sup>n</sup> hits are detected<br/>"
      + "&nbsp;&nbsp;&nbsp;&nbsp;plus at retention times that are predicted for species not identified by MS<sup>n</sup><html/>";
  public final static String RDI_GENERAL_ADDPOSITIONS = "<html>Are ther more chain positions possible than fatty acids<br/>"
      + "&nbsp;&nbsp;&nbsp;&nbsp;plus this is e.g. the case for DG class.<sup>n</sup><html/>";

  public final static String ISO_LABEL_LABELID = "<html>The designator for the label in the data</html>";
  public final static String ISO_LABEL_OMEGA = "<html>The omega position of the isotopically labeled chain</html>";
  public final static String ISO_LABEL_FORMULA = "<html>The difference in elemental composition caused by the label</html>";
  public final static String ISO_LABEL_RTSHIFT = "<html>The relative difference in retention time caused by the label</html>";
  public final static String ISO_DELETE = "<html>Deletes the label.</html>";
  public final static String ISO_ADD = "<html>Adds another row to enter information about a label</html>";
  
}
