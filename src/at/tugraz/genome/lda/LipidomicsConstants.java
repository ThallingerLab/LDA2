/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Objects;
import java.util.Properties;
import java.util.Vector;

import de.isas.mztab2.model.Contact;
import de.isas.mztab2.model.Instrument;
import de.isas.mztab2.model.Parameter;
import de.isas.mztab2.model.Publication;
import de.isas.mztab2.model.PublicationItem;
import de.isas.mztab2.model.PublicationItem.TypeEnum;
import de.isas.mztab2.model.Sample;
import de.isas.mztab2.model.SampleProcessing;
import at.tugraz.genome.lda.utils.Pair;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.xml.AbstractXMLSpectraReader;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class LipidomicsConstants
{
  public final static String EXCEL_MS_OH = "OH";
  public final static int EXCEL_NO_OH_INFO = -1;
  
  public final static String EXCEL_KEY = "Key";
  public final static String EXCEL_VALUE = "Value";
  
  public final static String EXCEL_MSN_SECTION_HEAD_FRAGMENTS = "Head group fragments";
  public final static String EXCEL_MSN_SECTION_HEAD_INTENSITIES = "Head group rules";
  public final static String EXCEL_MSN_SECTION_CHAIN_FRAGMENTS = "Chain fragments";
  public final static String EXCEL_MSN_SECTION_CHAIN_INTENSITIES = "Chain rules";
  public final static String EXCEL_MSN_SECTION_POSITION_INTENSITIES = "Position rules";
  
  public final static String EXCEL_MSN_FRAGMENT_NAME = "Name";
  public final static String EXCEL_MSN_FRAGMENT_OH = "OH";
  public final static String EXCEL_MSN_FRAGMENT_CHAIN_TYPE = "Chain type";  
  public final static String EXCEL_MSN_FRAGMENT_FORMULA = "Formula";
  public final static String EXCEL_MSN_FRAGMENT_CHARGE = "Charge";
  public final static String EXCEL_MSN_FRAGMENT_MSLEVEL = "msLevel";
  public final static String EXCEL_MSN_FRAGMENT_MZ = "Mz";
  public final static String EXCEL_MSN_FRAGMENT_MZ_TOLERANCE = "MzTolerance";
  public final static String EXCEL_MSN_FRAGMENT_AREA = "Area";
  public final static String EXCEL_MSN_FRAGMENT_PEAK = "Peak";
  public final static String EXCEL_MSN_FRAGMENT_TIME_LOWER = "LowerValley";
  public final static String EXCEL_MSN_FRAGMENT_TIME_UPPER = "UpperValley";
  public final static String EXCEL_MSN_FRAGMENT_MZ_LOWER = "LowMz";
  public final static String EXCEL_MSN_FRAGMENT_MZ_UPPER = "UpMz";
  public final static String EXCEL_MSN_FRAGMENT_ELLIPSE_TIME = "EllCentTime";
  public final static String EXCEL_MSN_FRAGMENT_ELLIPSE_MZ = "EllCentMz";
  public final static String EXCEL_MSN_FRAGMENT_ELLIPSE_TIME_RANGE = "EllStretchTime";
  public final static String EXCEL_MSN_FRAGMENT_ELLIPSE_MZ_RANGE = "EllStretchMz";

  public final static String EXCEL_MSN_INTENSITY_RULE = "Rule";
  public final static String EXCEL_MSN_INTENSITY_ORIGINAL = "Original";
  public final static String EXCEL_MSN_INTENSITY_VALUES = "Values";
  public final static String EXCEL_MSN_INTENSITY_MISSED = "Missed";
  
  public final static String EXCEL_HYDROXY_FA_PREFIX = "faOHEncoding_";
  public final static String EXCEL_HYDROXY_LCB_PREFIX = "lcbOHEncoding_";
  
  /** results shall be exported at the species level*/
  public final static short EXPORT_ANALYTE_TYPE_SPECIES = 0;
  /** results shall be exported at the chain level (no position)*/
  public final static short EXPORT_ANALYTE_TYPE_CHAIN = 1;
  /** results shall be exported at the chain position level*/
  public final static short EXPORT_ANALYTE_TYPE_POSITION = 2;
  
  /** to separate between the chains for naming of chain combinations*/
  public final static String CHAIN_COMBI_SEPARATOR = "<->";
  /** to indicate the number of OH sites after the type*/
  public final static String CHAIN_OH_INDEX_SEPARATOR = "^";
  /** to separate between the chain type and the name*/
  public final static String CHAIN_NAME_TYPE_SEPARATOR = "@";
  /** indicator that there will come a different type of linkage (alkyl/aleknyl)*/
  public final static String CHAIN_LINKAGE_INCLUSION_START = "{";
  /** indicator that there will end a different type of linkage (alkyl/aleknyl)*/
  public final static String CHAIN_LINKAGE_INCLUSION_STOP = "}";
  /** the prefix to indicate a fatty acyl/alkyl/alenyl chain*/
  public final static String CHAIN_TYPE_FA_NAME = "FA";
  /** the prefix to indicate an LCB chain*/
  public final static String CHAIN_TYPE_LCB_NAME = "LCB";
  /** the prefix to indicate that this is an fragment that was missed*/
  public final static short CHAIN_TYPE_MISSED = -2;
  /** the prefix to indicate a fatty acyl chain*/
  public final static short CHAIN_TYPE_NO_CHAIN = -1;
  /** the prefix to indicate a fatty acyl chain*/
  public final static short CHAIN_TYPE_FA_ACYL = 0;
  /** the prefix to indicate a fatty alkyl chain*/
  public final static short CHAIN_TYPE_FA_ALKYL = 1;
  /** the prefix to indicate a fatty alkyl chain*/
  public final static short CHAIN_TYPE_FA_ALKENYL = 2;  
  /** the prefix to indicate an LCB chain*/
  public final static short CHAIN_TYPE_LCB = 3;

  
  /** prefix for alkyl linked fatty acid chains */
  public final static String ALKYL_PREFIX = "O-";
  /** prefix for alkenyl linked fatty acid chains */
  public final static String ALKENYL_PREFIX = "P-";
  /** the separator to distinguish between several fatty acids*/
  public final static String CHAIN_SEPARATOR_NO_POS = "_";
  /** the separator to distinguish between several fatty acids, where the positions are known*/
  public final static String CHAIN_SEPARATOR_KNOWN_POS = "/";
  /** the separator between the number of carbon atoms and the double bonds*/  
  public final static String CHAIN_SEPARATOR_DBS = ":";
  /** separator used between human readable chain combinations to show that there are more than one possible position assignment*/
  public final static String CHAIN_COMBI_SEPARATOR_AMBIG_POS_OLD = ";";
  
  /** String starting the sn position assignment for proven positions*/
  public final static String SN_POSITION_START = "(sn-";
  /** String ending the sn position assignment for proven positions*/
  public final static String SN_POSITION_END = ")";
  
  /** String starting the omega position assignment*/
  public final static String OMEGA_POSITION_START = "(n-";
  /** String ending the omega position assignment*/
  public final static String OMEGA_POSITION_END = ")";
  
  /** String indicating that no FA is linked at this position*/
  public final static String NO_FA_LINKED = "-";
  
  /** the String for separating the OH numbers*/
  public final static String ALEX_OH_SEPARATOR = ";";
  /** prefix for alkylated chains*/
  //this was the old way of encoding - I am not sure whether this change will cause any errors
  ////public final static String ALEX_ALKYL_PREFIX = "O_";
  public final static String ALEX_ALKYL_PREFIX = "O-";
  /** prefix for alkenylated chains*/
//this was the old way of encoding - I am not sure whether this change will cause any errors
  //public final static String ALEX_ALKENYL_PREFIX = "P_";
  public final static String ALEX_ALKENYL_PREFIX = "P-";
  /** the String for separating the chains*/
  public final static String ALEX_CHAIN_SEPARATOR = "-";
  /** the String indicating an internal standard*/
  public final static String ALEX_IS_PREFIX = "IS ";

  
  
  /**chain-modification variables*/
  public final static String CHAIN_MOD_SEPARATOR =";";
  public final static String CHAIN_MOD_COLUMN_NAME ="PSM";
  
  public static String LDA_PROPERTIES_FILE = "LipidDataAnalyzer.properties";
  public static String MZTAB_PROPERTIES_FILE = "mzTab.properties";
  
  /** prefix for detectint addcut definitions in the MZTAB_PROPERTIES_FILE*/
  private final static String MZTAB_ADDUCT_PREFIX = "adduct_";
  
  private static LipidomicsConstants instance_ = null;
  private String ldaVersion_;
  private String rawFileName_;
    /** was this file quantified by using a Alex123 target list*/
  private boolean alexTargetlist_;
  /** lookup for classes whether MSn fragments were defined*/
  private Hashtable<String,Boolean> alexTargetlistUsed_;
  private String relativeMS1BasePeakCutoff_;
  private String currentMSMachine_;
  private float neutronMass_;
  private String basePeakDefaultCutoff_;
  private double massShift_;
  private int maxFileSizeForChromTranslationAtOnceInMB_;
  private float coarseChromMzTolerance_;
  private float chromSmoothRange_;
  private int chromSmoothRepeats_;
  private boolean isotopeCorrection_;
  private boolean removeFromOtherIsotopes_;
  private boolean respectIsotopicDistribution_;
  /** check if number of labels in species name corresponds with number of labels in chains*/
  private boolean checkChainLabelCombination_;
  /** minimum threshold in seconds to match predicted and measured retention times with high confidence */
  private int minimumThresholdForHighConfidenceRTMatch_;
  /** maximum threshold in seconds to match unambiguous predicted and measured retention times even if the match is not within the high confidence threshold */
  private int maximumThresholdForIntermediateConfidenceRTMatch_;
  private boolean useNoiseCutoff_;
  private float noiseCutoffDeviationValue_;
  private Float minimumRelativeIntensity_;
  private int scanStep_;
  private boolean use3D_;
  /** contains the acquired data only sparse MS1 data points*/
  private boolean sparseData_;
  private float profileMzRange_;
  private float profileTimeTolerance_;
  private float profileIntThreshold_;
  private float broaderProfileTimeTolerance_;   
  private float profileSmoothRange_;
  private int profileSmoothRepeats_;
  private int profileMeanSmoothRepeats_;
  private float profileMzMinRange_;
  private float profileSteepnessChange1_;
  private float profileSteepnessChange2_;
  private float profileIntensityCutoff1_;
  private float profileIntensityCutoff2_;
  private float profileGeneralIntCutoff_;
  private float profilePeakAcceptanceRange_;
  private float profileSmoothingCorrection_;
  private float profileMaxRange_;
  private float smallChromMzRange_;
  private int smallChromSmoothRepeats_;
  private int smallChromMeanSmoothRepeats_;
  private float smallChromSmoothRange_;
  private float smallChromIntensityCutoff_;
  private int broadChromSmoothRepeats_;
  private int broadChromMeanSmoothRepeats_;
  private int broadChromSmoothRange_;
  private float broadChromIntensityCutoff_;
  private float broadChromSteepnessChangeNoSmall_;
  private float broadChromIntensityCutoffNoSmall_;
  private float finalProbeTimeCompTolerance_;
  private float finalProbeMzCompTolerance_;
  private float overlapDistanceDeviationFactor_;
  private float overlapPossibleIntensityThreshold_;
  private float overlapSureIntensityThreshold_;
  private float overlapPeakDistanceDivisor_;
  private float overlapFullDistanceDivisor_;
  
  private int peakDiscardingAreaFactor_;
  private int isotopeInBetweenTime_;
  private float isoInBetweenAreaFactor_;
  private int isoInBetweenMaxTimeDistance_;
  private int isoNearNormalProbeTime_;
  private float relativeAreaCutoff_;
  private float relativeFarAreaCutoff_;
  private int relativeFarAreaTimeSpace_;
  private float relativeIsoInBetweenCutoff_;
  private int closePeakTimeTolerance_;
  private float twinInBetweenCutoff_;
  private float unionInBetweenCutoff_;

  private int chromMultiplicationFactorForInt_;
  private int chromLowestResolution_;
  
  private String threeDViewerDefaultTimeResolution_;
  private String threeDViewerDefaultMZResolution_;
  
  private boolean ms2_;
  private boolean useMostOverlappingIsotopeOnly_;
  
  private float ms2PrecursorTolerance_;
  private int ms2ChromMultiplicationFactorForInt_;
  private String threeDViewerMs2DefaultTimeResolution_;
  private String threeDViewerMs2DefaultMZResolution_;
  /** typically, for Waters files Mass++ is used; however, when you want to use msconvert, set this parameter to true*/
  private boolean useMsconvertForWaters_;
  
  /** the name of the instrument in PSI controlled vocabulary*/
  private Parameter mzTabInstrumentName_;
  /** the name of the ion source in PSI controlled vocabulary*/
  private Parameter mzTabInstrumentSource_;
  /** the name of the m/z analyzer in PSI controlled vocabulary*/
  private Parameter mzTabInstrumentAnalyzer_;
  /** the name of the instrument detector in PSI controlled vocabulary*/
  private Parameter mzTabInstrumentDetector_;
  
  private List<Contact> mzTabContacts_;
  /** for mztab-exporting: the name of the sample*/
  private Sample mzTabSample_;
  /** for mztab-exporting: the name sample processing steps*/
  private List<SampleProcessing> mzTabSampleProcessings_;
  /** for mztab-exporting: the list of publications the data refers to*/
  private List<Publication> mzTabPubs_;
  /** for mztab-exporting: the list of fragmentation methods the data refers to*/
  private List<Parameter> fragmethods_;
  /** for mztab-exporting: a lookup between the LDA adduct notation and the one proposed for the mzTab format*/
  private Hashtable<String,String> mzTabAdductLookup_;
  /** contains CV for SEP ontology*/
  private boolean cvSepRequired_;
  /** contains CV for CHMO ontology*/
  private boolean cvChmoRequired_;
  /** contains CV for NCBITaxon ontology*/
  private boolean cvNcbiTaxonRequired_;
  /** contains CV for CL ontology*/
  private boolean cvClRequired_;  
  /** contains CV for BTO ontology*/
  private boolean cvBtoRequired_;  
  /** contains CV for DOID ontology*/
  private boolean cvDoidRequired_;
  
  private boolean chromExportShowLegend_;
  
  //the following parameters are for MSn spectral analysis
  /** the +/- m/z tolerance for peak detection */
  private float ms2MzTolerance_;
  /** the +/- m/z tolerance for peak detection */  
  private int ms2MinIntsForNoiseRemoval_;
  /** relative cutoff threshold for fatty acid chain detection - it is in relation to the most intense chain */
  private double chainCutoffValue_;
  /** relative intensity cutoff for exclusion of isobars regarding spectrum coverage*/
  private float ms2IsobarSCExclusionRatio_;
  /** relative intensity cutoff for exclusion of isobars regarding spectrum coverage,
   * for areas that are farer away from a unique peak identification*/
  private float ms2IsobarSCFarExclusionRatio_;
  /** retention time in minutes that define a peak to be farer away, so that the ms2IsobarSCFarExclusionRatio can be used*/
  private float ms2IsobaricOtherRtDifference_;
  
  /** 0 when LC, 1 when shotgun, 2 when PRM*/
  private short shotgun_;
  /** the unit of the m/z values (available for shotgun only)*/
  private String mzUnit_;
  /** the way to acquire shotgun intensities (possible values are "mean", "median", "sum")*/
  private int shotgunProcessing_;
  /** should shotgun intensities be removed below a certain threshold*/
  private boolean shotgunIntensityRemoval_;
  /** the relative intensity cutoff for the removal*/
  private float shotgunRelIntCutoff_;
  /** the intermediate file format for the translation of RAW to chrom */
  private String intermediateFileFormat_;
  
  
  /** this is the conventional chromatography mode*/
  public final static short SHOTGUN_FALSE = 0;
  /** this is the shotgun mode*/
  public final static short SHOTGUN_TRUE = 1;
  /** this is the PRM mode*/
  public final static short SHOTGUN_PRM = 2;
  
  
  public final static String LDA_VERSION = "LDA-version";
  public final static String RAW_FILE = "rawFile";
  public final static String BASE_PEAK_CUTOFF = "basePeakCutoff";
  private final static String BASE_PEAK_CUTOFF_DEFAULT = "0.1";
  public final static String MASS_SHIFT = "massShift";
  private final static String MASS_SHIFT_DEFAULT = "0";
  private final static String MACHINE_NAME = "machineName";
  private final static String MACHINE_NAME_DEFAULT = "";
  private final static String NEUTRON_MASS = "neutronMass";
  private final static String NEUTRON_MASS_DEFAULT = "1.005";
  private final static String MAX_CHROM_AT_ONCE_MB = "maxFileSizeForChromTranslationAtOnce";
  private final static String MAX_CHROM_AT_ONCE_MB_DEFAULT = "50";
  private final static String COARSE_CHROM_MZ_TOL = "coarseChromMzTolerance";
  private final static String COARSE_CHROM_MZ_TOL_DEFAULT = "0.02";
  private final static String CHROM_SMOOTH_RANGE = "chromSmoothRange";
  private final static String CHROM_SMOOTH_RANGE_DEFAULT = "0.5";
  private final static String CHROM_SMOOTH_REPEATS = "chromSmoothRepeats";
  private final static String CHROM_SMOOTH_REPEATS_DEFAULT = "10";
  private final static String ISOTOPE_CORRECTION = "isotopeCorrection";
  private final static String ISOTOPE_CORRECTION_DEFAULT = "false";
  private final static String DEISOTOPE_BY_MOST_OVERLAPPING = "deisotopeByMostOverlappingOnly";
  private final static String DEISOTOPE_BY_MOST_OVERLAPPING_DEFAULT = "false";
  private final static String REMOVE_FROM_OTHER_ISOTOPE = "removeFromOtherIsotopes";
  private final static String REMOVE_FROM_OTHER_ISOTOPE_DEFAULT = "true";
  private final static String RESPECT_ISO_DISTRI = "respectIsotopicDistribution";
  private final static String RESPECT_ISO_DISTRI_DEFAULT = "true";
  private final static String CHECK_CHAIN_LABEL_COMBINATION = "checkChainLabelCombinationFromSpeciesName";
  private final static String CHECK_CHAIN_LABEL_COMBINATION_DEFAULT = "false";
  private final static String NOISE_CUTOFF = "useNoiseCutoff";
  private final static String NOISE_CUTOFF_DEFAULT = "false";
  private final static String NOISE_DEVIATION = "noiseCutoffDeviationValue";
  private final static String NOISE_DEVIATION_DEFAULT = "2";
  private final static String NOISE_MIN_INTENSITY = "minimumRelativeIntensity";
  private final static String SCAN_STEP = "scanStep";
  private final static String SCAN_STEP_DEFAULT = "1";
  private final static String USE_3D = "use3D";
  private final static String SPARSE_DATA = "sparseData";
  private final static String USE_3D_DEFAULT = "true";
  private final static String PROFILE_MZ_RANGE = "profileMzRangeExtraction";
  private final static String PROFILE_MZ_RANGE_DEFAULT = "0.2";
  private final static String PROFILE_TIME_TOL = "profileTimeTolerance";
  private final static String PROFILE_TIME_TOL_DEFAULT = "2";
  private final static String PROFILE_INT_THRESHOLD = "profileIntThreshold";
  private final static String PROFILE_INT_THRESHOLD_DEFAULT = "5";
  private final static String PROFILE_BROADER_TIME_TOL = "broaderProfileTimeTolerance";
  private final static String PROFILE_BROADER_TIME_TOL_DEFAULT = "3";
  private final static String PROFILE_SMOOTH_RANGE = "profileSmoothRange";
  private final static String PROFILE_SMOOTH_RANGE_DEFAULT = "0.02";
  private final static String PROFILE_SMOOTH_REPEATS = "profileSmoothRepeats";
  private final static String PROFILE_SMOOTH_REPEATS_DEFAULT = "5";
  private final static String PROFILE_MEAN_SMOOTH_REPEATS = "profileMeanSmoothRepeats";
  private final static String PROFILE_MEAN_SMOOTH_REPEATS_DEFAULT = "0";
  private final static String PROFILE_MZ_MIN_RANGE = "profileMzMinRange";
  private final static String PROFILE_MZ_MIN_RANGE_DEFAULT = "0.0";
  private final static String PROFILE_STEEPNESS_CHANGE_1 = "profileSteepnessChange1";
  private final static String PROFILE_STEEPNESS_CHANGE_1_DEFAULT = "1.5";
  private final static String PROFILE_STEEPNESS_CHANGE_2 = "profileSteepnessChange2";
  private final static String PROFILE_STEEPNESS_CHANGE_2_DEFAULT = "1.8";
  private final static String PROFILE_INT_CUTOFF_1 = "profileIntensityCutoff1";
  private final static String PROFILE_INT_CUTOFF_1_DEFAULT = "0.15";
  private final static String PROFILE_INT_CUTOFF_2 = "profileIntensityCutoff2";
  private final static String PROFILE_INT_CUTOFF_2_DEFAULT = "0.2";
  private final static String PROFILE_GENERAL_INT_CUTOFF = "profileGeneralIntCutoff";
  private final static String PROFILE_GENERAL_INT_CUTOFF_DEFAULT = "0.03";
  private final static String PROFILE_PEAK_ACCEPTANCE_RANGE = "profilePeakAcceptanceRange";
  private final static String PROFILE_PEAK_ACCEPTANCE_RANGE_DEFAULT = "0.034";
  private final static String PROFILE_SMOOTHING_CORRECTION = "profileSmoothingCorrection";
  private final static String PROFILE_SMOOTHING_CORRECTION_DEFAULT = "0.0";
  private final static String PROFILE_MZ_MAX_RANGE = "profileMaxRange";
  private final static String PROFILE_MZ_MAX_RANGE_DEFAULT = "0.02";
  private final static String SMALL_CHROM_MZ_RANGE = "smallChromMzRange";
  private final static String SMALL_CHROM_MZ_RANGE_DEFAULT = "0.002";
  private final static String SMALL_CHROM_SMOOTH_REPEATS = "smallChromSmoothRepeats";
  private final static String SMALL_CHROM_SMOOTH_REPEATS_DEFAULT = "15";
  private final static String SMALL_CHROM_MEAN_SMOOTH_REPEATS = "smallChromMeanSmoothRepeats";
  private final static String SMALL_CHROM_MEAN_SMOOTH_REPEATS_DEFAULT = "0";
  private final static String SMALL_CHROM_SMOOTH_RANGE = "smallChromSmoothRange";
  private final static String SMALL_CHROM_SMOOTH_RANGE_DEFAULT = "5";
  private final static String SMALL_CHROM_INT_CUTOFF = "smallChromIntensityCutoff";
  private final static String SMALL_CHROM_INT_CUTOFF_DEFAULT = "0.03";
  private final static String BROAD_CHROM_SMOOTH_REPEATS = "broadChromSmoothRepeats";
  private final static String BROAD_CHROM_SMOOTH_REPEATS_DEFAULT = "10";
  private final static String BROAD_CHROM_MEAN_SMOOTH_REPEATS = "broadChromMeanSmoothRepeats";
  private final static String BROAD_CHROM_MEAN_SMOOTH_REPEATS_DEFAULT = "0";
  private final static String BROAD_CHROM_SMOOTH_RANGE = "broadChromSmoothRange";
  private final static String BROAD_CHROM_SMOOTH_RANGE_DEFAULT = "5";
  private final static String BROAD_CHROM_INT_CUTOFF = "broadChromIntensityCutoff";
  private final static String BROAD_CHROM_INT_CUTOFF_DEFAULT = "0.0";
  private final static String BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL = "broadChromSteepnessChangeNoSmall";
  private final static String BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL_DEFAULT = "1.33";
  private final static String BROAD_CHROM_INT_CUTOFF_NO_SMALL = "broadChromIntensityCutoffNoSmall";
  private final static String BROAD_CHROM_INT_CUTOFF_NO_SMALL_DEFAULT = "0.05";
  private final static String FINAL_PROBE_TIME_COMP_TOL = "finalProbeTimeCompTolerance";
  private final static String FINAL_PROBE_TIME_COMP_TOL_DEFAULT = "0.1";
  private final static String FINAL_PROBE_MZ_COMP_TOL = "finalProbeMzCompTolerance";
  private final static String FINAL_PROBE_MZ_COMP_TOL_DEFAULT = "0.0005";
  private final static String OVERLAP_DIST_DEV_FACTOR = "overlapDistanceDeviationFactor";
  private final static String OVERLAP_DIST_DEV_FACTOR_DEFAULT = "1.5";
  private final static String OVERLAP_INT_THRESHOLD = "overlapPossibleIntensityThreshold";
  private final static String OVERLAP_INT_THRESHOLD_DEFAULT = "0.15";
  private final static String OVERLAP_INT_SURE_THRESHOLD = "overlapSureIntensityThreshold";
  private final static String OVERLAP_INT_SURE_THRESHOLD_DEFAULT = "0.7";
  private final static String OVERLAP_PEAK_DIST_DIVISOR = "overlapPeakDistanceDivisor";
  private final static String OVERLAP_PEAK_DIST_DIVISOR_DEFAULT = "3";
  private final static String OVERLAP_FULL_DIST_DIVISOR = "overlapFullDistanceDivisor";
  private final static String OVERLAP_FULL_DIST_DIVISOR_DEFAULT = "6";
  private final static String PEAK_DISCARD_AREA_FACTOR = "peakDiscardingAreaFactor";
  private final static String PEAK_DISCARD_AREA_FACTOR_DEFAULT = "1000";
  private final static String ISO_IN_BETWEEN_TIME = "isotopeInBetweenTime";
  private final static String ISO_IN_BETWEEN_TIME_DEFAULT = "30";
  private final static String ISO_IN_BETWEEN_AREA_FACTOR = "isoInBetweenAreaFactor";
  private final static String ISO_IN_BETWEEN_AREA_FACTOR_DEFAULT = "3";
  private final static String ISO_IN_BETWEEN_TIME_MAX = "isoInBetweenMaxTimeDistance";
  private final static String ISO_IN_BETWEEN_TIME_MAX_DEFAULT = "300";
  private final static String ISO_NEAR_NORMAL_PROBE_TIME = "isoNearNormalProbeTime";
  private final static String ISO_NEAR_NORMAL_PROBE_TIME_DEFAULT = "30";
  private final static String RELATIVE_AREA_CUTOFF = "relativeAreaCutoff";
  private final static String RELATIVE_AREA_CUTOFF_DEFAULT = "0.01";
  private final static String RELATIVE_AREA_FAR_CUTOFF = "relativeFarAreaCutoff";
  private final static String RELATIVE_AREA_FAR_CUTOFF_DEFAULT = "0.1";
  private final static String RELATIVE_AREA_FAR_TIME_SPACE = "relativeFarAreaTimeSpace";
  private final static String RELATIVE_AREA_FAR_TIME_SPACE_DEFAULT = "30";
  private final static String RELATIVE_ISO_INBETWEEN_CUTOFF = "relativeIsoInBetweenCutoff";
  private final static String RELATIVE_ISO_INBETWEEN_CUTOFF_DEFAULT = "0.5";
  private final static String PEAK_CLOSE_TIME_TOL = "closePeakTimeTolerance";
  private final static String PEAK_CLOSE_TIME_TOL_DEFAULT = "10";
  private final static String CHROM_MULT_FOR_INT = "chromMultiplicationFactorForInt";
  private final static String CHROM_MULT_FOR_INT_DEFAULT = "1000";
  private final static String CHROM_RESOLUTION_LOWEST = "chromLowestResolution";
  private final static String CHROM_RESOLUTION_LOWEST_DEFAULT = "1";
  private final static String VIEWER_TIME_RESOLUTION = "threeDViewerDefaultTimeResolution";
  private final static String VIEWER_TIME_RESOLUTION_DEFAULT = "5";
  private final static String VIEWER_MZ_RESOLUTION = "threeDViewerDefaultMZResolution";
  private final static String VIEWER_MZ_RESOLUTION_DEFAULT = "0.005";
  private final static String TWIN_INBETWEEN_CUTOFF = "twinInBetweenCutoff";
  private final static String TWIN_INBETWEEN_CUTOFF_DEFAULT = "0.95";
  private final static String UNION_INBETWEEN_CUTOFF = "unionInBetweenCutoff";
  private final static String UNION_INBETWEEN_CUTOFF_DEFAULT = "0.8";
  private final static String MS2 = "MS2";
  private final static String MS2_DEFAULT = "false";
  private final static String MS2_PRECURSOR_TOL = "ms2PrecursorTolerance";
  private final static String MS2_PRECURSOR_TOL_DEFAULT = "0.01";
  private final static String MS2_CHROM_MULT_FOR_INT = "ms2ChromMultiplicationFactorForInt";
  private final static String MS2_CHROM_MULT_FOR_INT_DEFAULT = "10";
  private final static String MS2_VIEWER_TIME_RESOLUTION = "threeDViewerMs2DefaultTimeResolution";
  private final static String MS2_VIEWER_TIME_RESOLUTION_DEFAULT = "1";
  private final static String MS2_VIEWER_MZ_RESOLUTION = "threeDViewerMs2DefaultMZResolution";
  private final static String MS2_VIEWER_MZ_RESOLUTION_DEFAULT = "1";
  private final static String MS2_MZ_TOL = "ms2MzTolerance";
  private final static String MS2_MIN_NOISE_REMOVAL = "ms2MinIntsForNoiseRemoval";
  private final static String MS2_ISOBAR_RATIO = "ms2IsobarSCExclusionRatio";
  private final static String MS2_ISOBAR_FAR_RATIO = "ms2IsobarSCFarExclusionRatio";
  private final static String MS2_ISOBAR_FAR_RT = "ms2IsobaricOtherRtDifference";
  private final static String MS2_MZ_TOL_DEFAULT = "0.2";
  private final static String MS2_MIN_NOISE_REMOVAL_DEFAULT = "100";
  private final static String CHAIN_CUTOFF = "chainCutoffValue";
  private final static String CHAIN_CUTOFF_DEFAULT = "0.01";
  private final static String ALEX_TARGETLIST = "alexTargetlist";
  private final static String USE_MSCONVERT_FOR_WATERS = "useMsconvertForWaters";
  private final static String MINIMUM_THRESHOLD_FOR_HIGH_CONFIDENCE_RT_MATCH = "minimumThresholdForHighConfidenceRTMatch";
  private final static String MINIMUM_THRESHOLD_FOR_HIGH_CONFIDENCE_RT_MATCH_DEFAULT = "4";
  private final static String MAXIMUM_THRESHOLD_FOR_INTERMEDIATE_CONFIDENCE_RT_MATCH = "maximumThresholdForIntermediateConfidenceRTMatch";
  private final static String MAXIMUM_THRESHOLD_FOR_INTERMEDIATE_CONFIDENCE_RT_MATCH_DEFAULT = "4";
  
  private final static String MZTAB_INSTRUMENT = "mzTabInstrumentName";
  private final static String MZTAB_IONSOURCE = "mzTabInstrumentIonsource";
  private final static String MZTAB_MSANALYZER = "mzTabInstrumentAnalyzer";
  private final static String MZTAB_DETECTOR = "mzTabInstrumentDetector";
  
  private final static String CHROM_EXPORT_LEGEND_IN_CHROM = "chromexportPaintLegendInChrom";
  
  /** the shotgun parameter - true/false/prm are allowed*/
  private final static String SHOTGUN = "shotgun";
  /** the input unit (for shotgan data only) - allowed are ppm, Da, Dalton, Th, Thompson, m/z*/
  private final static String MZUNIT = "mzUnit";
  /** the shotgun processing type - allowed are mean, median, sum*/
  private final static String SHOTGUN_PROCESSING = "shotgunProcessing";
  /** the way zeros are handled in calculating shotgun intensities ("false": zeros are used; "true": zeros are discarded; float value: relative cutoff to average without zeros)*/
  private final static String SHOTGUN_ZERO_HANDLING = "shotgunDiscardZeros";
  
  /** possible input parameter for MZUNIT*/ 
  public final static String MZUNIT_PPM = "ppm";
  /** possible input parameter for MZUNIT*/ 
  public final static String MZUNIT_DA = "Da";
  /** possible input parameter for MZUNIT*/
  private final static String MZUNIT_DALTON = "Dalton";
  /** possible input parameter for MZUNIT*/
  private final static String MZUNIT_THOMPSON = "Thompson";
  /** possible input parameter for MZUNIT*/
  private final static String MZUNIT_TH = "Th";
  /** possible input parameter for MZUNIT*/
  private final static String MZUNIT_MZ = "m/z";
  
  /** possible input parameter for SHOTGUN_PROCESSING*/
  private final static String SHOTGUN_PROCESSING_MEAN = "mean";
  /** possible input parameter for SHOTGUN_PROCESSING*/
  private final static String SHOTGUN_PROCESSING_MEDIAN = "median";
  /** possible input parameter for SHOTGUN_PROCESSING*/
  private final static String SHOTGUN_PROCESSING_SUM = "sum";
  /** the intermediate file format for the translation of RAW to chrom */
  private final static String INTERMEDIATE_FILE_FORMAT = "intermediateFileFormat";
  /** the default intermediate file format for the translation of RAW to chrom */
  private final static String INTERMEDIATE_FILE_FORMAT_DEFAULT = AbstractXMLSpectraReader.FILE_TYPE_MZ_XML;
  
  public static LipidomicsConstants getInstance() {
    if (instance_ == null) {
      instance_ = new LipidomicsConstants();
    }
    return instance_;
  }
  
  public LipidomicsConstants(boolean dummy){
    cvSepRequired_ = false;
    cvNcbiTaxonRequired_ = false;
    cvClRequired_ = false;
    cvBtoRequired_ = false;  
    cvDoidRequired_ = false;
    cvChmoRequired_ = false;
  }
  
  private LipidomicsConstants(){
    this(false);
    ldaVersion_ = Settings.VERSION;
    readConstantsFile(LDA_PROPERTIES_FILE);
    readmzTabConfFile(MZTAB_PROPERTIES_FILE);    
  }
  
  private void readConstantsFile(String fileName){
    try{
      File file = new File(fileName);
      System.out.println(file.getAbsolutePath());
      FileInputStream inNew = new FileInputStream(file);
      Properties properties = new Properties();
      properties.load(inNew);
      inNew.close();
      setVariables(properties);
    }catch(Exception e){
      e.printStackTrace();
    }    
  }
  
  public void setVariables(Properties properties) throws SettingsException{
    
    String cutoff = properties.getProperty(BASE_PEAK_CUTOFF, BASE_PEAK_CUTOFF_DEFAULT);
    try{
      Double.parseDouble(cutoff);
    }catch (NumberFormatException nfx){
      cutoff = BASE_PEAK_CUTOFF_DEFAULT;
    }
    massShift_ = 0d;
    String massShiftString = properties.getProperty(MASS_SHIFT,MASS_SHIFT_DEFAULT);
    try{
      massShift_ = Double.parseDouble(massShiftString)/1000d;
    } catch (NumberFormatException nfx){}

    currentMSMachine_ = properties.getProperty(MACHINE_NAME, MACHINE_NAME_DEFAULT);
    neutronMass_ = Float.parseFloat(properties.getProperty(NEUTRON_MASS,NEUTRON_MASS_DEFAULT));
    basePeakDefaultCutoff_ = cutoff;
    intermediateFileFormat_ = properties.getProperty(INTERMEDIATE_FILE_FORMAT,INTERMEDIATE_FILE_FORMAT_DEFAULT);
    maxFileSizeForChromTranslationAtOnceInMB_ = Integer.parseInt(properties.getProperty(MAX_CHROM_AT_ONCE_MB, MAX_CHROM_AT_ONCE_MB_DEFAULT));
    coarseChromMzTolerance_ = Float.parseFloat(properties.getProperty(COARSE_CHROM_MZ_TOL,COARSE_CHROM_MZ_TOL_DEFAULT));
    chromSmoothRange_ = Float.parseFloat(properties.getProperty(CHROM_SMOOTH_RANGE,CHROM_SMOOTH_RANGE_DEFAULT));
    chromSmoothRepeats_ = Integer.parseInt(properties.getProperty(CHROM_SMOOTH_REPEATS,CHROM_SMOOTH_REPEATS_DEFAULT));
    isotopeCorrection_ = false;
    String isotopeCorrectionString = properties.getProperty(ISOTOPE_CORRECTION,ISOTOPE_CORRECTION_DEFAULT);
    if (isotopeCorrectionString!=null && (isotopeCorrectionString.equalsIgnoreCase("yes")||isotopeCorrectionString.equalsIgnoreCase("true")))
      isotopeCorrection_ = true;
    useMostOverlappingIsotopeOnly_ = false;
    String mostOverlappingOnlyString = properties.getProperty(DEISOTOPE_BY_MOST_OVERLAPPING,DEISOTOPE_BY_MOST_OVERLAPPING_DEFAULT);
    if (mostOverlappingOnlyString!=null && (mostOverlappingOnlyString.equalsIgnoreCase("yes")||mostOverlappingOnlyString.equalsIgnoreCase("true")))
      useMostOverlappingIsotopeOnly_ = true;
    
    removeFromOtherIsotopes_ = true;
    String removeIsotopesString = properties.getProperty(REMOVE_FROM_OTHER_ISOTOPE,REMOVE_FROM_OTHER_ISOTOPE_DEFAULT);
    if (removeIsotopesString!=null && (removeIsotopesString.equalsIgnoreCase("no")||removeIsotopesString.equalsIgnoreCase("false")))
      removeFromOtherIsotopes_ = false;
    respectIsotopicDistribution_ = true;
    String respectIsotopesString = properties.getProperty(RESPECT_ISO_DISTRI,RESPECT_ISO_DISTRI_DEFAULT);
    if (respectIsotopesString!=null && (respectIsotopesString.equalsIgnoreCase("no")||respectIsotopesString.equalsIgnoreCase("false")))
      respectIsotopicDistribution_ = false;
    checkChainLabelCombination_ = false;
    String checkChainLabelCombinationString = properties.getProperty(CHECK_CHAIN_LABEL_COMBINATION,CHECK_CHAIN_LABEL_COMBINATION_DEFAULT);
    if (checkChainLabelCombinationString!=null && (checkChainLabelCombinationString.equalsIgnoreCase("yes")||checkChainLabelCombinationString.equalsIgnoreCase("true")))
      checkChainLabelCombination_ = true;
    minimumThresholdForHighConfidenceRTMatch_ = Integer.parseInt(properties.getProperty(MINIMUM_THRESHOLD_FOR_HIGH_CONFIDENCE_RT_MATCH,MINIMUM_THRESHOLD_FOR_HIGH_CONFIDENCE_RT_MATCH_DEFAULT));
    maximumThresholdForIntermediateConfidenceRTMatch_ = Integer.parseInt(properties.getProperty(MAXIMUM_THRESHOLD_FOR_INTERMEDIATE_CONFIDENCE_RT_MATCH,MAXIMUM_THRESHOLD_FOR_INTERMEDIATE_CONFIDENCE_RT_MATCH_DEFAULT));
    useNoiseCutoff_ = false;
    String cutoffString = properties.getProperty(NOISE_CUTOFF,NOISE_CUTOFF_DEFAULT);
    if (cutoffString!=null && (cutoffString.equalsIgnoreCase("yes")||cutoffString.equalsIgnoreCase("true")))
      useNoiseCutoff_ = true;
    noiseCutoffDeviationValue_ = Float.parseFloat(properties.getProperty(NOISE_DEVIATION,NOISE_DEVIATION_DEFAULT));
    minimumRelativeIntensity_ = null;
    String minIntString = properties.getProperty(NOISE_MIN_INTENSITY);
    if (minIntString!=null && minIntString.length()>0){
      try{
        minimumRelativeIntensity_ = new Float(minIntString);
      } catch (NumberFormatException nfx){}
    }
    scanStep_ = Integer.parseInt(properties.getProperty(SCAN_STEP,SCAN_STEP_DEFAULT));
    
    String use3DString = properties.getProperty(USE_3D,USE_3D_DEFAULT);
    use3D_ = false;
    if (use3DString!=null&&(use3DString.equalsIgnoreCase("true")||use3DString.equalsIgnoreCase("yes")))
      use3D_ = true;
    String sparseDataString = properties.getProperty(SPARSE_DATA,"false");
    if (sparseDataString!=null&&(sparseDataString.equalsIgnoreCase("true")||sparseDataString.equalsIgnoreCase("yes")))
      sparseData_ = true;
    profileMzRange_ = Float.parseFloat(properties.getProperty(PROFILE_MZ_RANGE,PROFILE_MZ_RANGE_DEFAULT))/2;
    profileTimeTolerance_ = Float.parseFloat(properties.getProperty(PROFILE_TIME_TOL,PROFILE_TIME_TOL_DEFAULT));
    profileIntThreshold_ = Float.parseFloat(properties.getProperty(PROFILE_INT_THRESHOLD,PROFILE_INT_THRESHOLD_DEFAULT));
    broaderProfileTimeTolerance_ = Float.parseFloat(properties.getProperty(PROFILE_BROADER_TIME_TOL,PROFILE_BROADER_TIME_TOL_DEFAULT));      
    profileSmoothRange_ = Float.parseFloat(properties.getProperty(PROFILE_SMOOTH_RANGE,PROFILE_SMOOTH_RANGE_DEFAULT));
    profileSmoothRepeats_ = Integer.parseInt(properties.getProperty(PROFILE_SMOOTH_REPEATS,PROFILE_SMOOTH_REPEATS_DEFAULT));
    profileMeanSmoothRepeats_ = Integer.parseInt(properties.getProperty(PROFILE_MEAN_SMOOTH_REPEATS,PROFILE_MEAN_SMOOTH_REPEATS_DEFAULT));
    profileMzMinRange_ = Float.parseFloat(properties.getProperty(PROFILE_MZ_MIN_RANGE,PROFILE_MZ_MIN_RANGE_DEFAULT));
    profileSteepnessChange1_ = Float.parseFloat(properties.getProperty(PROFILE_STEEPNESS_CHANGE_1,PROFILE_STEEPNESS_CHANGE_1_DEFAULT));
    profileSteepnessChange2_ = Float.parseFloat(properties.getProperty(PROFILE_STEEPNESS_CHANGE_2,PROFILE_STEEPNESS_CHANGE_2_DEFAULT));
    profileIntensityCutoff1_ = Float.parseFloat(properties.getProperty(PROFILE_INT_CUTOFF_1,PROFILE_INT_CUTOFF_1_DEFAULT));
    profileIntensityCutoff2_ = Float.parseFloat(properties.getProperty(PROFILE_INT_CUTOFF_2,PROFILE_INT_CUTOFF_2_DEFAULT));
    profileGeneralIntCutoff_ = Float.parseFloat(properties.getProperty(PROFILE_GENERAL_INT_CUTOFF,PROFILE_GENERAL_INT_CUTOFF_DEFAULT));
    profilePeakAcceptanceRange_ = Float.parseFloat(properties.getProperty(PROFILE_PEAK_ACCEPTANCE_RANGE,PROFILE_PEAK_ACCEPTANCE_RANGE_DEFAULT));
    profileSmoothingCorrection_ = Float.parseFloat(properties.getProperty(PROFILE_SMOOTHING_CORRECTION,PROFILE_SMOOTHING_CORRECTION_DEFAULT));
    profileMaxRange_ = Float.parseFloat(properties.getProperty(PROFILE_MZ_MAX_RANGE,PROFILE_MZ_MAX_RANGE_DEFAULT));
    smallChromMzRange_ = Float.parseFloat(properties.getProperty(SMALL_CHROM_MZ_RANGE,SMALL_CHROM_MZ_RANGE_DEFAULT));
    smallChromSmoothRepeats_ = Integer.parseInt(properties.getProperty(SMALL_CHROM_SMOOTH_REPEATS,SMALL_CHROM_SMOOTH_REPEATS_DEFAULT));
    smallChromMeanSmoothRepeats_ = Integer.parseInt(properties.getProperty(SMALL_CHROM_MEAN_SMOOTH_REPEATS,SMALL_CHROM_MEAN_SMOOTH_REPEATS_DEFAULT));
    smallChromSmoothRange_ = Float.parseFloat(properties.getProperty(SMALL_CHROM_SMOOTH_RANGE,SMALL_CHROM_SMOOTH_RANGE_DEFAULT));
    smallChromIntensityCutoff_ = Float.parseFloat(properties.getProperty(SMALL_CHROM_INT_CUTOFF,SMALL_CHROM_INT_CUTOFF_DEFAULT));
    broadChromSmoothRepeats_ = Integer.parseInt(properties.getProperty(BROAD_CHROM_SMOOTH_REPEATS,BROAD_CHROM_SMOOTH_REPEATS_DEFAULT));
    broadChromMeanSmoothRepeats_ = Integer.parseInt(properties.getProperty(BROAD_CHROM_MEAN_SMOOTH_REPEATS,BROAD_CHROM_MEAN_SMOOTH_REPEATS_DEFAULT));
    broadChromSmoothRange_ = Integer.parseInt(properties.getProperty(BROAD_CHROM_SMOOTH_RANGE,BROAD_CHROM_SMOOTH_RANGE_DEFAULT));
    broadChromIntensityCutoff_ = Float.parseFloat(properties.getProperty(BROAD_CHROM_INT_CUTOFF,BROAD_CHROM_INT_CUTOFF_DEFAULT));
    broadChromSteepnessChangeNoSmall_ = Float.parseFloat(properties.getProperty(BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL,BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL_DEFAULT));
    broadChromIntensityCutoffNoSmall_ = Float.parseFloat(properties.getProperty(BROAD_CHROM_INT_CUTOFF_NO_SMALL,BROAD_CHROM_INT_CUTOFF_NO_SMALL_DEFAULT));
    finalProbeTimeCompTolerance_ = Float.parseFloat(properties.getProperty(FINAL_PROBE_TIME_COMP_TOL,FINAL_PROBE_TIME_COMP_TOL_DEFAULT));
    finalProbeMzCompTolerance_ = Float.parseFloat(properties.getProperty(FINAL_PROBE_MZ_COMP_TOL,FINAL_PROBE_MZ_COMP_TOL_DEFAULT));
    overlapDistanceDeviationFactor_ = Float.parseFloat(properties.getProperty(OVERLAP_DIST_DEV_FACTOR,OVERLAP_DIST_DEV_FACTOR_DEFAULT));
    overlapPossibleIntensityThreshold_ = Float.parseFloat(properties.getProperty(OVERLAP_INT_THRESHOLD,OVERLAP_INT_THRESHOLD_DEFAULT));
    overlapSureIntensityThreshold_ = Float.parseFloat(properties.getProperty(OVERLAP_INT_SURE_THRESHOLD,OVERLAP_INT_SURE_THRESHOLD_DEFAULT));
    overlapPeakDistanceDivisor_ = Float.parseFloat(properties.getProperty(OVERLAP_PEAK_DIST_DIVISOR,OVERLAP_PEAK_DIST_DIVISOR_DEFAULT));
    overlapFullDistanceDivisor_ = Float.parseFloat(properties.getProperty(OVERLAP_FULL_DIST_DIVISOR,OVERLAP_FULL_DIST_DIVISOR_DEFAULT));
    peakDiscardingAreaFactor_= Integer.parseInt(properties.getProperty(PEAK_DISCARD_AREA_FACTOR,PEAK_DISCARD_AREA_FACTOR_DEFAULT));
    isotopeInBetweenTime_ = Integer.parseInt(properties.getProperty(ISO_IN_BETWEEN_TIME,ISO_IN_BETWEEN_TIME_DEFAULT));
    isoInBetweenAreaFactor_ = Float.parseFloat(properties.getProperty(ISO_IN_BETWEEN_AREA_FACTOR,ISO_IN_BETWEEN_AREA_FACTOR_DEFAULT));
    isoInBetweenMaxTimeDistance_ = Integer.parseInt(properties.getProperty(ISO_IN_BETWEEN_TIME_MAX,ISO_IN_BETWEEN_TIME_MAX_DEFAULT));
    isoNearNormalProbeTime_ = Integer.parseInt(properties.getProperty(ISO_NEAR_NORMAL_PROBE_TIME,ISO_NEAR_NORMAL_PROBE_TIME_DEFAULT));
    relativeAreaCutoff_ = Float.parseFloat(properties.getProperty(RELATIVE_AREA_CUTOFF,RELATIVE_AREA_CUTOFF_DEFAULT));
    relativeFarAreaCutoff_ = Float.parseFloat(properties.getProperty(RELATIVE_AREA_FAR_CUTOFF,RELATIVE_AREA_FAR_CUTOFF_DEFAULT));
    relativeFarAreaTimeSpace_ = Integer.parseInt(properties.getProperty(RELATIVE_AREA_FAR_TIME_SPACE,RELATIVE_AREA_FAR_TIME_SPACE_DEFAULT));
    relativeIsoInBetweenCutoff_ = Float.parseFloat(properties.getProperty(RELATIVE_ISO_INBETWEEN_CUTOFF,RELATIVE_ISO_INBETWEEN_CUTOFF_DEFAULT));
    closePeakTimeTolerance_ = Integer.parseInt(properties.getProperty(PEAK_CLOSE_TIME_TOL,PEAK_CLOSE_TIME_TOL_DEFAULT));
    chromMultiplicationFactorForInt_ = Integer.parseInt(properties.getProperty(CHROM_MULT_FOR_INT,CHROM_MULT_FOR_INT_DEFAULT));
    chromLowestResolution_ = Integer.parseInt(properties.getProperty(CHROM_RESOLUTION_LOWEST,CHROM_RESOLUTION_LOWEST_DEFAULT));
    threeDViewerDefaultTimeResolution_ = properties.getProperty(VIEWER_TIME_RESOLUTION,VIEWER_TIME_RESOLUTION_DEFAULT);
    try{
      Double.parseDouble(threeDViewerDefaultTimeResolution_);
    } catch (NumberFormatException nfx){threeDViewerDefaultTimeResolution_ = VIEWER_TIME_RESOLUTION_DEFAULT;}
    threeDViewerDefaultMZResolution_ = properties.getProperty(VIEWER_MZ_RESOLUTION,VIEWER_MZ_RESOLUTION_DEFAULT);
    try{
      Double.parseDouble(threeDViewerDefaultMZResolution_);
    } catch (NumberFormatException nfx){threeDViewerDefaultTimeResolution_ = VIEWER_MZ_RESOLUTION_DEFAULT;}
    twinInBetweenCutoff_ = Float.parseFloat(properties.getProperty(TWIN_INBETWEEN_CUTOFF,TWIN_INBETWEEN_CUTOFF_DEFAULT));
    unionInBetweenCutoff_ = Float.parseFloat(properties.getProperty(UNION_INBETWEEN_CUTOFF,UNION_INBETWEEN_CUTOFF_DEFAULT));
    
    String ms2String = properties.getProperty(MS2,MS2_DEFAULT);
    ms2_ = false;
    if (ms2String!=null&&(ms2String.equalsIgnoreCase("true")||ms2String.equalsIgnoreCase("yes")))
      ms2_ = true;
    ms2PrecursorTolerance_ = Float.parseFloat(properties.getProperty(MS2_PRECURSOR_TOL,MS2_PRECURSOR_TOL_DEFAULT));
    ms2ChromMultiplicationFactorForInt_ = Integer.parseInt(properties.getProperty(MS2_CHROM_MULT_FOR_INT,MS2_CHROM_MULT_FOR_INT_DEFAULT));
    threeDViewerMs2DefaultTimeResolution_ = properties.getProperty(MS2_VIEWER_TIME_RESOLUTION ,MS2_VIEWER_TIME_RESOLUTION_DEFAULT);
    try{
      Float.parseFloat(threeDViewerMs2DefaultTimeResolution_);
    } catch (NumberFormatException nfx){threeDViewerMs2DefaultTimeResolution_ = MS2_VIEWER_TIME_RESOLUTION_DEFAULT;}
    threeDViewerMs2DefaultMZResolution_ = properties.getProperty(MS2_VIEWER_MZ_RESOLUTION,MS2_VIEWER_MZ_RESOLUTION_DEFAULT);
    try{
      Float.parseFloat(threeDViewerMs2DefaultMZResolution_);
    } catch (NumberFormatException nfx){threeDViewerMs2DefaultMZResolution_ = MS2_VIEWER_MZ_RESOLUTION_DEFAULT;}
    ms2MzTolerance_ = Float.parseFloat(properties.getProperty(MS2_MZ_TOL,MS2_MZ_TOL_DEFAULT));
    try {
      ms2MinIntsForNoiseRemoval_ = Integer.parseInt(properties.getProperty(MS2_MIN_NOISE_REMOVAL,MS2_MIN_NOISE_REMOVAL_DEFAULT));
    } catch (NumberFormatException nfx) {
      throw new SettingsException("Invalid input for "+MS2_MIN_NOISE_REMOVAL+"! The input must be in integer format!");
    }
    
    ms2IsobarSCExclusionRatio_ = 0f;
    String ms2IsobarSCExclusionRatioString = properties.getProperty(MS2_ISOBAR_RATIO);
    if (ms2IsobarSCExclusionRatioString!=null && ms2IsobarSCExclusionRatioString.length()>0){
      try{
        ms2IsobarSCExclusionRatio_ = Float.parseFloat(ms2IsobarSCExclusionRatioString);
      } catch (NumberFormatException nfx){}
    }
    ms2IsobarSCFarExclusionRatio_ = 0f;
    String ms2IsobarSCFarExclusionRatioString = properties.getProperty(MS2_ISOBAR_FAR_RATIO);
    if (ms2IsobarSCFarExclusionRatioString!=null && ms2IsobarSCFarExclusionRatioString.length()>0){
      try{
        ms2IsobarSCFarExclusionRatio_ = Float.parseFloat(ms2IsobarSCFarExclusionRatioString);
      } catch (NumberFormatException nfx){}
    }
    ms2IsobaricOtherRtDifference_ = Float.MAX_VALUE;
    String ms2IsobaricOtherRtDifferenceString = properties.getProperty(MS2_ISOBAR_FAR_RT);
    if (ms2IsobaricOtherRtDifferenceString!=null && ms2IsobaricOtherRtDifferenceString.length()>0){
      try{
        ms2IsobaricOtherRtDifference_ = Float.parseFloat(ms2IsobaricOtherRtDifferenceString);
      } catch (NumberFormatException nfx){}
    }
        
    chainCutoffValue_ = Double.parseDouble(properties.getProperty(CHAIN_CUTOFF,CHAIN_CUTOFF_DEFAULT));
    
    
    String chromExportString = properties.getProperty(CHROM_EXPORT_LEGEND_IN_CHROM);
    chromExportShowLegend_ = false;
    if (chromExportString!=null&&(chromExportString.equalsIgnoreCase("true")||chromExportString.equalsIgnoreCase("yes")))
      chromExportShowLegend_ = true;
    
    String alexTargetlistString = properties.getProperty(ALEX_TARGETLIST,"false");
    alexTargetlist_ = false;
    if (alexTargetlistString!=null&&(alexTargetlistString.equalsIgnoreCase("true")||alexTargetlistString.equalsIgnoreCase("yes")))
      alexTargetlist_ = true;
    alexTargetlistUsed_ = new Hashtable<String,Boolean>();
    
    String useMsconvertForWatersString = properties.getProperty(USE_MSCONVERT_FOR_WATERS,"false");
    useMsconvertForWaters_ = false;
    if (useMsconvertForWatersString!=null&&(useMsconvertForWatersString.equalsIgnoreCase("true")||useMsconvertForWatersString.equalsIgnoreCase("yes")))
      useMsconvertForWaters_ = true;
    
    mzTabInstrumentName_ = extractEBIParam(MZTAB_INSTRUMENT,properties);
    mzTabInstrumentSource_ = extractEBIParam(MZTAB_IONSOURCE,properties);
    mzTabInstrumentAnalyzer_ = extractEBIParam(MZTAB_MSANALYZER,properties);
    mzTabInstrumentDetector_ = extractEBIParam(MZTAB_DETECTOR,properties);
    
    String shotgunString = properties.getProperty(SHOTGUN,"false");
    shotgun_ = SHOTGUN_FALSE;
    if (shotgunString!=null&&(shotgunString.equalsIgnoreCase("true")||shotgunString.equalsIgnoreCase("yes")))
      shotgun_ = SHOTGUN_TRUE;
    if (shotgunString!=null&&(shotgunString.equalsIgnoreCase("prm")))
      shotgun_ = SHOTGUN_PRM;
    
    String mzUnitString = properties.getProperty(MZUNIT,MZUNIT_PPM);
    mzUnit_= MZUNIT_PPM;
    if (shotgun_==SHOTGUN_FALSE || mzUnitString.equalsIgnoreCase(MZUNIT_DA) || mzUnitString.equalsIgnoreCase(MZUNIT_DALTON) ||
        mzUnitString.equalsIgnoreCase(MZUNIT_THOMPSON) || mzUnitString.equalsIgnoreCase(MZUNIT_TH) ||
        mzUnitString.equalsIgnoreCase(MZUNIT_MZ))
      mzUnit_ = MZUNIT_DA;
    shotgunProcessing_ = LipidomicsAnalyzer.SHOTGUN_TYPE_MEAN;
    if (shotgun_==SHOTGUN_TRUE){
      String processingString = properties.getProperty(SHOTGUN_PROCESSING);
      if (processingString==null)
        throw new SettingsException("When there is shotgun data, the \""+SHOTGUN_PROCESSING+"\" has to be set!");
      else if (processingString.equalsIgnoreCase(SHOTGUN_PROCESSING_MEAN))
        shotgunProcessing_ = LipidomicsAnalyzer.SHOTGUN_TYPE_MEAN;
      else if (processingString.equalsIgnoreCase(SHOTGUN_PROCESSING_MEDIAN))
        shotgunProcessing_ = LipidomicsAnalyzer.SHOTGUN_TYPE_MEDIAN;
      else if (processingString.equalsIgnoreCase(SHOTGUN_PROCESSING_SUM))
        shotgunProcessing_ = LipidomicsAnalyzer.SHOTGUN_TYPE_SUM;
      else
        throw new SettingsException("The shotgun processing type \""+processingString+"\" is unknown! Please use "+SHOTGUN_PROCESSING_MEAN+", "+SHOTGUN_PROCESSING_MEDIAN+", or "+SHOTGUN_PROCESSING_SUM+" instead!");
    }
    
    String shotgunZeroHandlingString = properties.getProperty(SHOTGUN_ZERO_HANDLING,"false");
    this.shotgunIntensityRemoval_ = false;
    this.shotgunRelIntCutoff_ = -1f;
    if (shotgunZeroHandlingString.equalsIgnoreCase("false") || shotgunZeroHandlingString.equalsIgnoreCase("no"))
      ;
    else if (shotgunZeroHandlingString.equalsIgnoreCase("true") || shotgunZeroHandlingString.equalsIgnoreCase("yes")) {
      this.shotgunIntensityRemoval_ = true;
      this.shotgunRelIntCutoff_ = 0f;
    } else {
      try {
        shotgunRelIntCutoff_ = (float)StaticUtils.readPercentPermilleValue(shotgunZeroHandlingString);
        shotgunIntensityRemoval_ = true;
      } catch (NumberFormatException nfx) {
        throw new SettingsException("Invalid input for "+SHOTGUN_ZERO_HANDLING+"! The following values are allowed: true, false, or any float format that might be followed by % or \u2030!");
      }
    }
      
  }
  
  /**
   * @return the intermediate file format for the translation of RAW to chrom
   */
  public static String getIntermediateFileFormat()
  {
    getInstance();
    return instance_.intermediateFileFormat_;
  }
  
  /**
   * @return the file pieces for the file translation to the chrom file
   */
  public static int getmMaxFileSizeForChromTranslationAtOnceInMB()
  {
    getInstance();
    return instance_.maxFileSizeForChromTranslationAtOnceInMB_;
  }
  
  /**
   * @return the m/z tolerance for the first coarse chromatogram
   */
  public static float getCoarseChromMzTolerance(float mz)
  {
    getInstance();
    return getCorrectMzTolerance(instance_.coarseChromMzTolerance_,instance_.mzUnit_,mz);
  }
  
  public static double getMassShift()
  {
    getInstance();
    return instance_.massShift_;
  }
  
  public static float getChromSmoothRange()
  {
    getInstance();
    return instance_.chromSmoothRange_;
  }

  public static int getChromSmoothRepeats()
  {
    getInstance();
    return instance_.chromSmoothRepeats_;
  }

  public static boolean isotopicCorrection()
  {
    getInstance();
    return instance_.isotopeCorrection_;
  }
  
  public static boolean removeIfOtherIsotopePresent()
  {
    getInstance();
    return instance_.removeFromOtherIsotopes_;
  }

  public static boolean removeIfDistriDoesNotFit()
  {
    getInstance();
    return instance_.respectIsotopicDistribution_;
  }
  
  public static boolean useNoiseCutoff()
  {
    getInstance();
    return instance_.useNoiseCutoff_;
  }
  
  /**
   * 
   * @return check if number of labels in species name corresponds with number of labels in chains
   */
  public static boolean checkChainLabelCombination()
  {
    getInstance();
    return instance_.checkChainLabelCombination_;
  }
  
  public static int getMinimumThresholdForHighConfidenceRTMatch() 
  {
    getInstance();
    return instance_.minimumThresholdForHighConfidenceRTMatch_;
  }
  
  public static int getMaximumThresholdForIntermediateConfidenceRTMatch() 
  {
    getInstance();
    return instance_.maximumThresholdForIntermediateConfidenceRTMatch_;
  }
  
  public static float getNoiseCutoffDeviationValue()
  {
    getInstance();
    return instance_.noiseCutoffDeviationValue_;
  }
  
  public static Float getMinimumRelativeIntensity()
  {
    getInstance();
    return instance_.minimumRelativeIntensity_;
  }
  

  public static int getScanStep()
  {
    getInstance();
    return instance_.scanStep_;
  }
  
  /**
   * @return if 3D-method should be used for automated quantification
   */
  public static boolean use3D()
  {
    getInstance();
    return instance_.use3D_;
  }
  
  /**
   * @return the m/z tolerance for the first coarse chromatogram
   */
  public static float getProfileMzRange()
  {
    getInstance();
    return instance_.profileMzRange_;
  }
  
  /**
   * @return the time tolerance of the profile
   */
  public static float getProfileTimeTolerance_()
  {
    getInstance();
    return instance_.profileTimeTolerance_;
  }
  
  /**
   * @return the intensity threshold for the extraction of a new the profile
   */
  public static float getProfileIntThreshold_()
  {
    getInstance();
    return instance_.profileIntThreshold_;
  }
  
  /**
   * @return the range for the broader profile extraction
   */
  public static float getBroaderProfileTimeTolerance_()
  {
    getInstance();
    return instance_.broaderProfileTimeTolerance_;
  }
  
  /**
   * @return the m/z smoothing range for the profile
   */
  public static float getProfileSmoothRange()
  {
    getInstance();
    return instance_.profileSmoothRange_;
  }
  
  /**
   * @return the amount of smoothing repeats for the profile
   */
  public static int getProfileSmoothRepeats()
  {
    getInstance();
    return instance_.profileSmoothRepeats_;
  }
  /**
   * @return the amount of smoothing repeats for the profile
   */
  public static int getProfileMeanSmoothRepeats()
  {
    getInstance();
    return instance_.profileMeanSmoothRepeats_;
  }  
  /**
   * @return the amount of smoothing repeats for the profile
   */
  public static float getProfileMzMinRange()
  {
    getInstance();
    return instance_.profileMzMinRange_;
  }
  /**
   * @return steepness change1 threshold for profile
   */
  public static float getProfileSteepnessChange1()
  {
    getInstance();
    return instance_.profileSteepnessChange1_;
  }
  
  /**
   * @return steepness change2 threshold for profile
   */
  public static float getProfileSteepnessChange2()
  {
    getInstance();
    return instance_.profileSteepnessChange2_;
  }  
  
  /**
   * @return steepness intensity cutoff 1 threshold for profile
   */
  public static float getProfileIntensityCutoff1()
  {
    getInstance();
    return instance_.profileIntensityCutoff1_;
  }  

  /**
   * @return steepness intensity cutoff 2 threshold for profile
   */
  public static float getProfileIntensityCutoff2()
  {
    getInstance();
    return instance_.profileIntensityCutoff2_;
  }  

  /**
   * @return steepness general intensity cutoff threshold for profile
   */
  public static float getProfileGeneralIntCutoff()
  {
    getInstance();
    return instance_.profileGeneralIntCutoff_;
  }
  
  /**
   * @return steepness general intensity cutoff threshold for profile
   */
  public static float getProfilePeakAcceptanceRange()
  {
    getInstance();
    return instance_.profilePeakAcceptanceRange_;
  }
  
  /**
   * @return m/z correction for profile smoothing
   */
  public static float getProfileSmoothingCorrection()
  {
    getInstance();
    return instance_.profileSmoothingCorrection_;
  }
  
  
  /**
   * @return the maximum size of one profile
   */
  public static float getProfileMaxRange()
  {
    getInstance();
    return instance_.profileMaxRange_;
  }
  
  /**
   * @return m/z range for the small chromatogram
   */
  public static float getSmallChromMzRange()
  {
    getInstance();
    return instance_.smallChromMzRange_;
  } 
  
  /**
   * @return smoothing repeats for the small chromatogram
   */
  public static int getSmallChromSmoothRepeats()
  {
    getInstance();
    return instance_.smallChromSmoothRepeats_;
  }
  
  /**
   * @return smoothing repeats for the small chromatogram
   */
  public static int getSmallChromMeanSmoothRepeats()
  {
    getInstance();
    return instance_.smallChromMeanSmoothRepeats_;
  }
  
  /**
   * @return smoothing repeats for the small chromatogram
   */
  public static float getSmallChromSmoothRange()
  {
    getInstance();
    return instance_.smallChromSmoothRange_;
  } 
  
  /**
   * @return intensity-cutoff for the small chromatogram
   */
  public static float getSmallChromIntensityCutoff()
  {
    getInstance();
    return instance_.smallChromIntensityCutoff_;
  }
  
  /**
   * @return smoothing repeats for the broad chromatogram
   */
  public static int getBroadChromSmoothRepeats()
  {
    getInstance();
    return instance_.broadChromSmoothRepeats_;
  }
  
  /**
   * @return smoothing repeats for the broad chromatogram
   */
  public static int getBroadChromMeanSmoothRepeats()
  {
    getInstance();
    return instance_.broadChromMeanSmoothRepeats_;
  }
  
  /**
   * @return smoothing range for the broad chromatogram
   */
  public static float getBroadChromSmoothRange()
  {
    getInstance();
    return instance_.broadChromSmoothRange_;
  } 
  
  /**
   * @return intensity-cutoff for the broad chromatogram
   */
  public static float getBroadChromIntensityCutoff()
  {
    getInstance();
    return instance_.broadChromIntensityCutoff_;
  }  
  
  /**
   * @return steepness-change threshold for the broad chromatogram if one of the 
   * small borders is not valid
   */
  public static float getBroadChromSteepnessChangeNoSmall()
  {
    getInstance();
    return instance_.broadChromSteepnessChangeNoSmall_;
  }
  
  /**
   * @return intensity-cutoff for the broad chromatogram if one of the 
   * small borders is not valid
   */
  public static float getBroadIntensityCutoffNoSmall()
  {
    getInstance();
    return instance_.broadChromIntensityCutoffNoSmall_;
  }
  
  /**
   * @return the tolerance or a comparison of the final probe in time direction
   */
  public static float getFinalProbeTimeCompTolerance()
  {
    getInstance();
    return instance_.finalProbeTimeCompTolerance_;
  }
  
  /**
   * @return the tolerance for a comparison of the final probe in m/z direction
   */
  public static float getFinalProbeMzCompTolerance()
  {
    getInstance();
    return instance_.finalProbeMzCompTolerance_;
  }

  
  /**
   * @return the distance tolerance factor for the detection of an overlap;
   * the comparison is done between the original probe and the newly calculated one
   */
  public static float getOverlapDistanceDeviationFactor()
  {
    getInstance();
    return instance_.overlapDistanceDeviationFactor_;
  }
  
  /**
   * @return a multiplication factor for the comparison of intensities
   * determines if an overlap is possible
   */
  public static float getOverlapPossibleIntensityThreshold()
  {
    getInstance();
    return instance_.overlapPossibleIntensityThreshold_;
  }

  /**
   * @return a multiplication factor for the comparison of intensities
   * determines if an overlap is sure
   */
  public static float getOverlapSureIntensityThreshold()
  {
    getInstance();
    return instance_.overlapSureIntensityThreshold_;
  }
  /**
   * @return a divisor to check if the peak of an isotopic peak is near
   * the other peak to assume an overlap
   * this distance is between peak and one border
   */
  public static float getOverlapPeakDistanceDivisor()
  {
    getInstance();
    return instance_.overlapPeakDistanceDivisor_;
  }
  /**
   * @return a divisor to check if the peak of an isotopic peak is near
   * the other peak to assume an overlap
   * this distance is between the borders
   */
  public static float getOverlapFullDistanceDivisor()
  {
    getInstance();
    return instance_.overlapFullDistanceDivisor_;
  }
  
  /**
   * @return peaks that are relatively smaller to the highest found peak than this factor are discarded
   */
  public static int getPeakDiscardingAreaFactor()
  {
    getInstance();
    return instance_.peakDiscardingAreaFactor_;
  }
  /**
   * @return allowed time distance in seconds to check if there is an isotope between the two peaks
   */
  public static int getIsotopeInBetweenTime()
  {
    getInstance();
    return instance_.isotopeInBetweenTime_;
  }
  /**
   * @return if an isotope is in between the smaller one of the two peaks is regarded to belong to be a fragment
   * of the other isotope if the area of the smaller one is isoInBetweenAreaFactor times smaller than
   * the bigger one
   */
  public static float getIsoInBetweenAreaFactor()
  {
    getInstance();
    return instance_.isoInBetweenAreaFactor_;
  }
  
  public static int getIsoInBetweenMaxTimeDistance()
  {
    getInstance();
    return instance_.isoInBetweenMaxTimeDistance_;
  }
  
  /**
   * @return this is the time distance in seconds where an isotopic peak is still regarded as interesting
   */
  public static int getIsoNearNormalProbeTime()
  {
    getInstance();
    return instance_.isoNearNormalProbeTime_;
  }
  /**
   * @return after the isotopic peaks are thrown out a higher threshold is applied to discard small peaks 
   * (relative to the highest found one)
   */
  public static float getRelativeAreaCutoff()
  {
    getInstance();
    return instance_.relativeAreaCutoff_;
  }
  /**
   * @return intensities that are close to the highest one can be members of a twin peak, hence they are
   * not discarded so quickly
   * this is a relative threshold to discard peaks that are farer away
   */
  public static float getRelativeFarAreaCutoff()
  {
    getInstance();
    return instance_.relativeFarAreaCutoff_;
  }
  /**
   * @return intensities that are close to the highest one can be members of a twin peak, hence they are
   * not discarded so quickly
   * this is the time distance in seconds to define if a peak is near the highest one
   */
  public static int getRelativeFarAreaTimeSpace()
  {
    getInstance();
    return instance_.relativeFarAreaTimeSpace_;
  }
  public static float getTwinInBetweenCutoff()
  {
    getInstance();
    return instance_.twinInBetweenCutoff_;
  }
  public static float getUnionInBetweenCutoff()
  {
    getInstance();
    return instance_.unionInBetweenCutoff_;
  }
  /**
   * @return if still there is an isotopic peak in between the smaller peak is normally discarded
   * if the ratio between smallerPeak/higherPeak the smaller one is not discarded since
   * it cannot be clearly identified which one of the peaks is the correct one
   */
  public static float getRelativeIsoInBetweenCutoff()
  {
    getInstance();
    return instance_.relativeIsoInBetweenCutoff_;
  }
  
  /**
   * @return this is the time distance in seconds defines if in a last step close peaks should be united
   */
  public static int getClosePeakTimeTolerance()
  {
    getInstance();
    return instance_.closePeakTimeTolerance_;
  }
  
  public static String getBasePeakDefaultCutoff()
  {
    getInstance();
    return instance_.basePeakDefaultCutoff_;
  }
  public static int getChromMultiplicationFactorForInt(){
    getInstance();
    return instance_.chromMultiplicationFactorForInt_;
  }

  public static int getChromLowestResolution(){
    getInstance();
    return instance_.chromLowestResolution_;
  }

  public static String getThreeDViewerDefaultTimeResolution(){
    getInstance();
    return instance_.threeDViewerDefaultTimeResolution_;
  }
  
  public static String getThreeDViewerDefaultMZResolution(){
    getInstance();
    return instance_.threeDViewerDefaultMZResolution_;
  }
  
  /**
   * @return if MS2 should be used for analyisis
   */
  public static boolean isMS2()
  {
    getInstance();
    return instance_.ms2_;
  }
  
  public static float getMs2PrecursorTolerance(){
    getInstance();
    return instance_.ms2PrecursorTolerance_;
  }
  
  public static int getMs2ChromMultiplicationFactorForInt(){
    getInstance();
    return instance_.ms2ChromMultiplicationFactorForInt_;
  }
  
  public static String getThreeDViewerMs2DefaultTimeResolution_(){
    getInstance();
    return instance_.threeDViewerMs2DefaultTimeResolution_;
  }

  public static String getThreeDViewerMs2DefaultMZResolution(){
    getInstance();
    return instance_.threeDViewerMs2DefaultMZResolution_;
  }
  
  public static boolean isChromExportShowLegend()
  {
    getInstance();
    return instance_.chromExportShowLegend_;
  }
  
  /**
   * 
   * @return the +/- m/z tolerance for peak detection
   */
  public static float getMs2MzTolerance(){
    getInstance();
    return instance_.ms2MzTolerance_;
  }
  
  /**
   * 
   * @return minimum number of detected signals to start a noise removal
   */
  public static int getMs2MinIntsForNoiseRemoval(){
    getInstance();
    return instance_.ms2MinIntsForNoiseRemoval_;
  }
  
  /**
   * 
   * @return relative intensity cutoff for exclusion of isobars regarding spectrum coverage
   */
  public static float getMs2IsobarSCExclusionRatio(){
    getInstance();
    return instance_.ms2IsobarSCExclusionRatio_;
  }

  
  /**
   * 
   * @return relative intensity cutoff for exclusion of isobars regarding spectrum coverage, for areas that are farer away from a unique peak identification
   */
  public static float getMs2IsobarSCFarExclusionRatio(){
    getInstance();
    return instance_.ms2IsobarSCFarExclusionRatio_;
  }

  /**
   * 
   * @return retention time in minutes that define a peak to be farer away, so that the ms2IsobarSCFarExclusionRatio can be used
   */
  public static float getMs2IsobaricOtherRtDifference(){
    getInstance();
    return instance_.ms2IsobaricOtherRtDifference_;
  }
  
  /**
   * 
   * @return 0 when LC, 1 when shotgun, 2 when PRM
   */
  public static short isShotgun(){
    getInstance();
    return instance_.shotgun_;
  }
  
  /**
   * 
   * @return the unit type of the m/z value
   */
  public static String getMzUnit(){
    getInstance();
    return instance_.mzUnit_;    
  }
  
  /**
   * 
   * @return the shotgun processing type according to the definitions in LipidomicsAnalyzer
   */
  public static int getShogunProcessing(){
    getInstance();
    return instance_.shotgunProcessing_;    
  }
  
  /**
   * 
   * @return should shotgun intensities be removed below a certain threshold
   */
  public static boolean isShotgunIntensityRemoval(){
    getInstance();
    return instance_.shotgunIntensityRemoval_;    
  }

  /**
   * 
   * @return the relative intensity cutoff for the removal
   */
  public static float getShotgunRelIntCutoff(){
    getInstance();
    return instance_.shotgunRelIntCutoff_;    
  }

  /**
   * 
   * @return relative cutoff threshold for fatty acid chain detection - it is in relation to the most intense chain
   */
  public static double getChainCutoffValue(){
    getInstance();
    return instance_.chainCutoffValue_;
  }
  
  public static String getCurrentMSMachine(){
    getInstance();
    return instance_.currentMSMachine_;
  }
  
  public static float getNeutronMass(){
    getInstance();
    return instance_.neutronMass_;
  }
  
  public static boolean useMostOverlappingIsotopeOnly(){
    getInstance();    
    return instance_.useMostOverlappingIsotopeOnly_;
  }
  
  /**
   * 
   * @return contains the acquired data only sparse MS1 data points
   */
  public static boolean isSparseData(){
    getInstance();
    return instance_.sparseData_;
  }
  
  public static Instrument getMzTabInstrument(){
    getInstance();
    if (instance_.mzTabInstrumentName_!=null || instance_.mzTabInstrumentSource_!=null || instance_.mzTabInstrumentAnalyzer_!=null || instance_.mzTabInstrumentDetector_!=null){
      Instrument instrument = new Instrument().id(1);
      instrument.addAnalyzerItem(instance_.mzTabInstrumentAnalyzer_);
      instrument.setDetector(instance_.mzTabInstrumentDetector_);
      instrument.setName(instance_.mzTabInstrumentName_);
      instrument.setSource(instance_.mzTabInstrumentSource_);
      return instrument;
    }
    return null;
  }
  
  public static List<Contact> getMzTabContacts(){
    getInstance();
    return instance_.mzTabContacts_;
  }

  /**
   * 
   * @return for mzTab-export: list of sample processing steps defined in the MZTAB_PROPERTIES_FILE
   */
  public static List<SampleProcessing> getMzTabSampleprocessings(){
    getInstance();
    return instance_.mzTabSampleProcessings_;
  }

  /**
   * 
   * @return for mzTab-export: list of publications this data refers to
   */
  public static List<Publication> getMzTabPublications(){
    getInstance();
    return instance_.mzTabPubs_;
  }
  
  /**
   * 
   * @return for mzTab-export: list of fragmentation methods used for the data
   */
  public static List<Parameter> getFragmentationMethods(){
    getInstance();
    return instance_.fragmethods_;
  }

  /**
   * 
   * @return for mzTab-export: the sample this data originates of
   */
  public static Sample getMzTabSample(){
    getInstance();
    return instance_.mzTabSample_;    
  }
  
  /**
   * 
   * @return true when this is a shotgun instance
   */
  public short getShotgun(){
    getInstance();
    return this.shotgun_;
  }
  
  /**
   * checks the adduct lookup and returns the corresponding mzTab notation when present, null otherwise
   * @param ldaAdduct the LDA adduct notation
   * @return the mzTab notation
   */
  public static String getMzTabAdduct(String ldaAdduct){
    getInstance();
    String mzTabAdduct = null;
    if (instance_.mzTabAdductLookup_.containsKey(ldaAdduct))
      mzTabAdduct = instance_.mzTabAdductLookup_.get(ldaAdduct);
    return mzTabAdduct;
  }
  
  /**
   * 
   * @return contains CV for SEP ontology
   */
  public static boolean isCvSepRequired(){
    getInstance();
    return instance_.cvSepRequired_;    
  }

  /**
   * 
   * @return contains CV for NCBITaxon ontology
   */
  public static boolean isCvNcbiTaxonRequired(){
    getInstance();
    return instance_.cvNcbiTaxonRequired_;
  }

  /**
   * 
   * @return contains CV for CL ontology
   */
  public static boolean isClRequired(){
    getInstance();
    return instance_.cvClRequired_;
  }
  
  /**
   * 
   * @return contains CV for BTO ontology
   */
  public static boolean isCvBtoRequired(){
    getInstance();
    return instance_.cvBtoRequired_;
  }

  /**
   * 
   * @return contains CV for DOID ontology
   */
  public static boolean isCvDoidRequired(){
    getInstance();
    return instance_.cvDoidRequired_;
  }
  
  /**
   * 
   * @return contains CV for CHMO ontology
   */
  public static boolean isCvChmoRequired(){
    getInstance();
    return instance_.cvChmoRequired_;    
  }
  
  public static void switchToOtherConfFile(File newConfFile){
    getInstance();
    instance_.readConstantsFile(newConfFile.getAbsolutePath());
  }
    
  public static void switchToOtherDefaultConfFile(File newConfFile) throws IOException{
    getInstance();
    int chunkSize = 1024;
    InputStream in = new BufferedInputStream(new FileInputStream(newConfFile));
    OutputStream out = new BufferedOutputStream(new FileOutputStream(LDA_PROPERTIES_FILE));
    byte[] buffer = new byte[chunkSize];
    int len;
    while ((len = in.read(buffer)) > 0) {
      out.write(buffer, 0, len);
    }
    in.close();
    out.close();
    instance_.readConstantsFile(LDA_PROPERTIES_FILE);
  }
  
  /**
   * parses a properties entry and returns the corresponding parameter in PSI notation
   * @param key the key of the property
   * @param properties the properties read from the file
   * @return the corresponding parameter in PSI notation
   */
  private Parameter extractEBIParam(String key, Properties properties){
    Parameter param = null;
    String paramString = properties.getProperty(key,""); 
    if (paramString!=null && paramString.length()>0 && paramString.indexOf(",")!=-1){
      Vector<String> split = splitByComma(paramString);
      if (split.size()==4){
        //TODO: I used her "" instead of null, since the jmzTab lib has a small bug
        String cvLabel = "";
        if (split.get(0)!=null && split.get(0).length()>0){
          cvLabel = split.get(0);
          if (!cvSepRequired_ && cvLabel.equalsIgnoreCase("SEP")) cvSepRequired_ = true;
          else if (!cvNcbiTaxonRequired_ && cvLabel.equalsIgnoreCase("NCBITaxon")) cvNcbiTaxonRequired_ = true;
          else if (!cvClRequired_ && cvLabel.equalsIgnoreCase("CL")) cvClRequired_ = true;
          else if (!cvBtoRequired_ && cvLabel.equalsIgnoreCase("BTO")) cvBtoRequired_ = true;
          else if (!cvDoidRequired_ && cvLabel.equalsIgnoreCase("DOID")) cvDoidRequired_ = true;
          if (!cvChmoRequired_ && cvLabel.equalsIgnoreCase("CHMO")) cvChmoRequired_ = true;
        }
        //TODO: I used her "" instead of null, since the jmzTab lib has a small bug
        String accession = "";
        if (split.get(1)!=null && split.get(1).length()>0) accession = split.get(1);
        //TODO: I used her "" instead of null, since the jmzTab lib has a small bug
        String name = "";
        if (split.get(2)!=null && split.get(2).length()>0) name = split.get(2);
        String value = null;
        if (split.get(3)!=null && split.get(3).length()>0) value = split.get(3);
        param = new Parameter().cvLabel(cvLabel).cvAccession(accession).name(name).value(value);
      }
    }
    return param;
  }
  
  private static Vector<String> splitByComma(String str){
    Vector<String> splitString = new  Vector<String>();
    String clone = new String(str);
    int pos = -1;
    while ((pos = clone.indexOf(","))!=-1){
      splitString.add(clone.substring(0,pos).trim());
      clone = clone.substring(pos+1);
    }
    splitString.add(clone.trim());
    return splitString;
  }
  
  private static int countComma(String str){
    String clone = new String(str);
    int pos = -1;
    int count = 0;
    while ((pos = clone.indexOf(","))!=-1){
      clone = clone.substring(pos+1);
      count++;
    }
    return count;
  }

  
  private void readmzTabConfFile(String fileName){
    try{
      File file = new File(fileName);
      mzTabAdductLookup_ = new Hashtable<String,String>();
      if (file.exists()){
        FileInputStream inNew = new FileInputStream(file);
        Properties properties = new Properties();
        properties.load(inNew);
        inNew.close();
        mzTabContacts_ = new ArrayList<Contact>();
        String propBase = "contact_";
        String propBase2 = "";
        int count = 1;
        String contactString = "";
        while ((contactString = properties.getProperty(propBase+String.valueOf(count),""))!=null && contactString.length()>0){
          count++;
          if (countComma(contactString)<2) continue;
          String name = contactString.substring(0,contactString.indexOf(","));
          contactString = contactString.substring(contactString.indexOf(",")+1);
          String email = contactString.substring(0,contactString.indexOf(","));
          String affiliation = contactString.substring(contactString.indexOf(",")+1);
          Contact contact = new Contact();
          contact.setId(count-1);
          contact.setAffiliation(affiliation);
          contact.setEmail(email);
          contact.setName(name);
          mzTabContacts_.add(contact);
        }
        mzTabSampleProcessings_ = new ArrayList<SampleProcessing>();
        String processingString = "";
        propBase = "sampleprocessing_";
        count = 1;
        while ((processingString = properties.getProperty(propBase+String.valueOf(count),""))!=null && processingString.length()>0){
          count++;
          int commas = countComma(processingString);
          if (countComma(processingString)<1) continue;
          SampleProcessing processing = new SampleProcessing().id(count-1);
          Parameter procDetails = null;
          if (commas==1){
            procDetails = new Parameter().name(processingString.substring(0,processingString.indexOf(","))).value(processingString.substring(processingString.indexOf(",")+1));
          }else{
            procDetails = extractEBIParam(propBase+String.valueOf(count-1),properties);
          }
          processing.addSampleProcessingItem(procDetails);
          mzTabSampleProcessings_.add(processing);
        }
        mzTabPubs_ = new ArrayList<Publication>();
        String publicationString1 = "";
        String publicationString2 = "";
        propBase = "pubmedid_";
        propBase2 = "doi_";
        count = 1;
        while (((publicationString1 = properties.getProperty(propBase+String.valueOf(count),""))!=null && publicationString1.length()>0) |
               ((publicationString2 = properties.getProperty(propBase2+String.valueOf(count),""))!=null && publicationString2.length()>0)){
          count++;
          
          Publication pub = new Publication().id(count-1);
          if (publicationString1!=null && publicationString1.length()>0){
            PublicationItem item = new PublicationItem().accession(publicationString1).type(TypeEnum.PUBMED);
            pub.addPublicationItemsItem(item);
          }
          if (publicationString2!=null && publicationString2.length()>0){
            PublicationItem item = new PublicationItem().accession(publicationString2).type(TypeEnum.DOI);
            pub.addPublicationItemsItem(item);            
          }
          mzTabPubs_.add(pub);
        }
        mzTabSample_ = null;
        Parameter mzTabSpecies = extractEBIParam("species",properties);
        Parameter mzTabTissue = extractEBIParam("tissue",properties);
        Parameter mzTabCelltype = extractEBIParam("celltype",properties);
        List<Parameter> speciesDiseases = extractParameterList("species_disease_",properties);
        if (mzTabSpecies!=null || mzTabTissue!=null || mzTabCelltype!=null){
          mzTabSample_ = new Sample().id(1);
          if (mzTabSpecies!=null){
            mzTabSample_.setSpecies(new Vector<Parameter>());
            mzTabSample_.getSpecies().add(mzTabSpecies);
          }
          if (mzTabTissue!=null){
            mzTabSample_.setTissue((new Vector<Parameter>()));
            mzTabSample_.getTissue().add(mzTabTissue);
          }
          if (mzTabCelltype!=null){
            mzTabSample_.setCellType((new Vector<Parameter>()));
            mzTabSample_.getCellType().add(mzTabCelltype);
          }
          if (speciesDiseases.size()>0)
            mzTabSample_.setDisease(speciesDiseases);
        }
        fragmethods_ = extractParameterList("fragmentationmethod_",properties);
        
        String key;
        String value;
        String adductKey;
        for (Object keyObject : properties.keySet()){
          if (!(keyObject instanceof String))
            continue;
          key = (String) keyObject;
          if (key.startsWith(MZTAB_ADDUCT_PREFIX)){
            adductKey = key.substring(MZTAB_ADDUCT_PREFIX.length());
            value = properties.getProperty(key);
            mzTabAdductLookup_.put(adductKey,value);
          }
        }
      }
    }catch(Exception e){
      e.printStackTrace();
    }
  }
  
  /**
   * extracts a sorted list of parameters from the properties which start with the propBase followed by a consecutive number starting with one  
   * @param propBase the base of the properties to extract
   * @param properties available properties
   * @return sorted list of extracted parameters
   */
  private List<Parameter> extractParameterList(String propBase, Properties properties){
    List<Parameter> parameters = new ArrayList<Parameter>();
    int count = 1;
    String value = "";
    while ((value = properties.getProperty(propBase+String.valueOf(count),""))!=null && value.length()>0){    
      int commas = countComma(value);
      if (countComma(value)<1) continue;
      Parameter param = null;
      if (commas==1){
        param = new Parameter().name(value.substring(0,value.indexOf(","))).value(value.substring(value.indexOf(",")+1));
      }else{
        param = extractEBIParam(propBase+String.valueOf(count),properties);
      }
      parameters.add(param);
      count++;
    }
    return parameters;
  }
  
  
  public List<Pair<String,String>> getPropertyRowList(HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) {
    List<Pair<String,String>> propertyRows = new ArrayList<Pair<String,String>>();
    propertyRows.add(new Pair<String,String>(EXCEL_KEY,EXCEL_VALUE));
    propertyRows.add(new Pair<String,String>(LDA_VERSION,ldaVersion_));
    if (rawFileName_!=null){
      propertyRows.add(new Pair<String,String>(RAW_FILE,rawFileName_));
    }
    propertyRows.add(new Pair<String,String>(MACHINE_NAME,currentMSMachine_));
    propertyRows.add(new Pair<String,String>(NEUTRON_MASS,String.valueOf(neutronMass_)));
    propertyRows.add(new Pair<String,String>(COARSE_CHROM_MZ_TOL,String.valueOf(coarseChromMzTolerance_)));
    propertyRows.add(new Pair<String,String>(MS2,String.valueOf(ms2_)));
    if (relativeMS1BasePeakCutoff_!=null){
      propertyRows.add(new Pair<String,String>(BASE_PEAK_CUTOFF,String.valueOf(relativeMS1BasePeakCutoff_)));
    }
    propertyRows.add(new Pair<String,String>(MASS_SHIFT,String.valueOf(massShift_)));
    propertyRows.add(new Pair<String,String>(VIEWER_TIME_RESOLUTION,String.valueOf(threeDViewerDefaultTimeResolution_)));
    propertyRows.add(new Pair<String,String>(VIEWER_MZ_RESOLUTION,String.valueOf(threeDViewerDefaultMZResolution_)));
    propertyRows.add(new Pair<String,String>(MS2_PRECURSOR_TOL,String.valueOf(ms2PrecursorTolerance_)));
    propertyRows.add(new Pair<String,String>(MS2_MZ_TOL,String.valueOf(ms2MzTolerance_)));
    propertyRows.add(new Pair<String,String>(MS2_MIN_NOISE_REMOVAL,String.valueOf(ms2MinIntsForNoiseRemoval_)));
    propertyRows.add(new Pair<String,String>(MS2_ISOBAR_RATIO,String.valueOf(ms2IsobarSCExclusionRatio_)));
    propertyRows.add(new Pair<String,String>(MS2_ISOBAR_FAR_RATIO,String.valueOf(ms2IsobarSCFarExclusionRatio_)));
    propertyRows.add(new Pair<String,String>(MS2_ISOBAR_FAR_RT,String.valueOf(ms2IsobaricOtherRtDifference_)));
    propertyRows.add(new Pair<String,String>(CHAIN_CUTOFF,String.valueOf(chainCutoffValue_)));
    propertyRows.add(new Pair<String,String>(MS2_CHROM_MULT_FOR_INT,String.valueOf(ms2ChromMultiplicationFactorForInt_)));
    propertyRows.add(new Pair<String,String>(MS2_VIEWER_TIME_RESOLUTION,threeDViewerMs2DefaultTimeResolution_));
    propertyRows.add(new Pair<String,String>(MS2_VIEWER_MZ_RESOLUTION,threeDViewerMs2DefaultMZResolution_));
    propertyRows.add(new Pair<String,String>(MAX_CHROM_AT_ONCE_MB,String.valueOf(maxFileSizeForChromTranslationAtOnceInMB_)));
    propertyRows.add(new Pair<String,String>(CHROM_MULT_FOR_INT,String.valueOf(chromMultiplicationFactorForInt_)));
    propertyRows.add(new Pair<String,String>(CHROM_RESOLUTION_LOWEST,String.valueOf(chromLowestResolution_)));
    if (shotgun_>SHOTGUN_FALSE){
      propertyRows.add(new Pair<String,String>(SHOTGUN,String.valueOf(shotgun_)));
      propertyRows.add(new Pair<String,String>(MZUNIT,mzUnit_));
    }
    if (shotgun_==SHOTGUN_TRUE){
      String shotgunProcessingString = null;
      switch(shotgunProcessing_) {
        case LipidomicsAnalyzer.SHOTGUN_TYPE_MEAN:
          shotgunProcessingString = SHOTGUN_PROCESSING_MEAN; break;
        case LipidomicsAnalyzer.SHOTGUN_TYPE_MEDIAN:
          shotgunProcessingString = SHOTGUN_PROCESSING_MEDIAN; break;
        case LipidomicsAnalyzer.SHOTGUN_TYPE_SUM:
          shotgunProcessingString = SHOTGUN_PROCESSING_SUM; break;
      }
      propertyRows.add(new Pair<String,String>(SHOTGUN_PROCESSING,shotgunProcessingString));
      
      String shotgunZeroHandlingString = "false";
      if (shotgunRelIntCutoff_>0)
        shotgunZeroHandlingString = String.valueOf(shotgunRelIntCutoff_);
      else if (shotgunRelIntCutoff_==0)
        shotgunZeroHandlingString = "true";
      propertyRows.add(new Pair<String,String>(SHOTGUN_ZERO_HANDLING,shotgunZeroHandlingString));
    } else {
      propertyRows.add(new Pair<String,String>(CHROM_SMOOTH_RANGE,String.valueOf(chromSmoothRange_)));
      propertyRows.add(new Pair<String,String>(CHROM_SMOOTH_REPEATS,String.valueOf(chromSmoothRepeats_)));
    }
    if (shotgun_==SHOTGUN_FALSE){
      propertyRows.add(new Pair<String,String>(USE_3D,String.valueOf(use3D_)));
      propertyRows.add(new Pair<String,String>(ISOTOPE_CORRECTION,String.valueOf(isotopeCorrection_)));
      propertyRows.add(new Pair<String,String>(REMOVE_FROM_OTHER_ISOTOPE,String.valueOf(removeFromOtherIsotopes_)));
      propertyRows.add(new Pair<String,String>(RESPECT_ISO_DISTRI,String.valueOf(respectIsotopicDistribution_)));
      propertyRows.add(new Pair<String,String>(CHECK_CHAIN_LABEL_COMBINATION,String.valueOf(checkChainLabelCombination_)));
      propertyRows.add(new Pair<String,String>(NOISE_CUTOFF,String.valueOf(useNoiseCutoff_)));
      propertyRows.add(new Pair<String,String>(NOISE_DEVIATION,String.valueOf(noiseCutoffDeviationValue_)));
      if (minimumRelativeIntensity_!=null) {
        propertyRows.add(new Pair<String,String>(NOISE_MIN_INTENSITY,minimumRelativeIntensity_.toString()));
      }
      propertyRows.add(new Pair<String,String>(SCAN_STEP,String.valueOf(scanStep_)));
      propertyRows.add(new Pair<String,String>(PROFILE_MZ_RANGE,String.valueOf(profileMzRange_*2)));
      propertyRows.add(new Pair<String,String>(PROFILE_TIME_TOL,String.valueOf(profileTimeTolerance_)));
      propertyRows.add(new Pair<String,String>(PROFILE_INT_THRESHOLD,String.valueOf(profileIntThreshold_)));
      propertyRows.add(new Pair<String,String>(PROFILE_BROADER_TIME_TOL,String.valueOf(broaderProfileTimeTolerance_)));
      propertyRows.add(new Pair<String,String>(PROFILE_SMOOTH_RANGE,String.valueOf(profileSmoothRange_)));
      propertyRows.add(new Pair<String,String>(PROFILE_SMOOTH_REPEATS,String.valueOf(profileSmoothRepeats_)));
      propertyRows.add(new Pair<String,String>(PROFILE_MEAN_SMOOTH_REPEATS,String.valueOf(profileMeanSmoothRepeats_)));
      propertyRows.add(new Pair<String,String>(PROFILE_MZ_MIN_RANGE,String.valueOf(profileMzMinRange_)));
      propertyRows.add(new Pair<String,String>(PROFILE_STEEPNESS_CHANGE_1,String.valueOf(profileSteepnessChange1_)));
      propertyRows.add(new Pair<String,String>(PROFILE_STEEPNESS_CHANGE_2,String.valueOf(profileSteepnessChange2_)));
      propertyRows.add(new Pair<String,String>(PROFILE_INT_CUTOFF_1,String.valueOf(profileIntensityCutoff1_)));
      propertyRows.add(new Pair<String,String>(PROFILE_INT_CUTOFF_2,String.valueOf(profileIntensityCutoff2_)));
      propertyRows.add(new Pair<String,String>(PROFILE_GENERAL_INT_CUTOFF,String.valueOf(profileGeneralIntCutoff_)));
      propertyRows.add(new Pair<String,String>(PROFILE_PEAK_ACCEPTANCE_RANGE,String.valueOf(profilePeakAcceptanceRange_)));
      propertyRows.add(new Pair<String,String>(PROFILE_SMOOTHING_CORRECTION,String.valueOf(profileSmoothingCorrection_)));
      propertyRows.add(new Pair<String,String>(PROFILE_MZ_MAX_RANGE,String.valueOf(profileMaxRange_)));
      propertyRows.add(new Pair<String,String>(SMALL_CHROM_MZ_RANGE,String.valueOf(smallChromMzRange_)));
      propertyRows.add(new Pair<String,String>(SMALL_CHROM_SMOOTH_REPEATS,String.valueOf(smallChromSmoothRepeats_)));
      propertyRows.add(new Pair<String,String>(SMALL_CHROM_MEAN_SMOOTH_REPEATS,String.valueOf(smallChromMeanSmoothRepeats_)));
      propertyRows.add(new Pair<String,String>(SMALL_CHROM_SMOOTH_RANGE,String.valueOf(smallChromSmoothRange_)));
      propertyRows.add(new Pair<String,String>(SMALL_CHROM_INT_CUTOFF,String.valueOf(smallChromIntensityCutoff_)));
      propertyRows.add(new Pair<String,String>(BROAD_CHROM_SMOOTH_REPEATS,String.valueOf(broadChromSmoothRepeats_)));
      propertyRows.add(new Pair<String,String>(BROAD_CHROM_MEAN_SMOOTH_REPEATS,String.valueOf(broadChromMeanSmoothRepeats_)));
      propertyRows.add(new Pair<String,String>(BROAD_CHROM_SMOOTH_RANGE,String.valueOf(broadChromSmoothRange_)));
      propertyRows.add(new Pair<String,String>(BROAD_CHROM_INT_CUTOFF,String.valueOf(broadChromIntensityCutoff_)));
      propertyRows.add(new Pair<String,String>(BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL,String.valueOf(broadChromSteepnessChangeNoSmall_)));
      propertyRows.add(new Pair<String,String>(FINAL_PROBE_TIME_COMP_TOL,String.valueOf(finalProbeTimeCompTolerance_)));
      propertyRows.add(new Pair<String,String>(FINAL_PROBE_MZ_COMP_TOL,String.valueOf(finalProbeMzCompTolerance_)));
      propertyRows.add(new Pair<String,String>(OVERLAP_DIST_DEV_FACTOR,String.valueOf(overlapDistanceDeviationFactor_)));
      propertyRows.add(new Pair<String,String>(OVERLAP_INT_THRESHOLD,String.valueOf(overlapPossibleIntensityThreshold_)));
      propertyRows.add(new Pair<String,String>(OVERLAP_INT_SURE_THRESHOLD,String.valueOf(overlapSureIntensityThreshold_)));
      propertyRows.add(new Pair<String,String>(OVERLAP_PEAK_DIST_DIVISOR,String.valueOf(overlapPeakDistanceDivisor_)));
      propertyRows.add(new Pair<String,String>(OVERLAP_FULL_DIST_DIVISOR,String.valueOf(overlapFullDistanceDivisor_)));
    }
    if (shotgun_!=SHOTGUN_TRUE){
      propertyRows.add(new Pair<String,String>(PEAK_DISCARD_AREA_FACTOR,String.valueOf(peakDiscardingAreaFactor_)));
    }
    if (shotgun_==SHOTGUN_FALSE){
      propertyRows.add(new Pair<String,String>(ISO_IN_BETWEEN_TIME,String.valueOf(isotopeInBetweenTime_)));
      propertyRows.add(new Pair<String,String>(ISO_IN_BETWEEN_AREA_FACTOR,String.valueOf(isoInBetweenAreaFactor_)));
      propertyRows.add(new Pair<String,String>(ISO_NEAR_NORMAL_PROBE_TIME,String.valueOf(isoNearNormalProbeTime_)));
    }
    if (shotgun_!=SHOTGUN_TRUE){
      propertyRows.add(new Pair<String,String>(RELATIVE_AREA_CUTOFF,String.valueOf(relativeAreaCutoff_)));
      propertyRows.add(new Pair<String,String>(RELATIVE_AREA_FAR_CUTOFF,String.valueOf(relativeFarAreaCutoff_)));
      propertyRows.add(new Pair<String,String>(RELATIVE_AREA_FAR_TIME_SPACE,String.valueOf(relativeFarAreaTimeSpace_)));
    }
    if (shotgun_==SHOTGUN_FALSE){
      propertyRows.add(new Pair<String,String>(RELATIVE_ISO_INBETWEEN_CUTOFF,String.valueOf(relativeIsoInBetweenCutoff_)));
      propertyRows.add(new Pair<String,String>(ISO_IN_BETWEEN_TIME_MAX,String.valueOf(isoInBetweenMaxTimeDistance_)));
      propertyRows.add(new Pair<String,String>(PEAK_CLOSE_TIME_TOL,String.valueOf(closePeakTimeTolerance_)));
      propertyRows.add(new Pair<String,String>(TWIN_INBETWEEN_CUTOFF,String.valueOf(twinInBetweenCutoff_)));
      propertyRows.add(new Pair<String,String>(UNION_INBETWEEN_CUTOFF,String.valueOf(unionInBetweenCutoff_)));
      propertyRows.add(new Pair<String,String>(SPARSE_DATA,String.valueOf(sparseData_)));
    }
    if (mzTabInstrumentName_!=null){
      propertyRows.add(new Pair<String,String>(MZTAB_INSTRUMENT,getStringFromMzTabParameter(mzTabInstrumentName_)));
    }
    if (mzTabInstrumentSource_!=null){
      propertyRows.add(new Pair<String,String>(MZTAB_IONSOURCE,getStringFromMzTabParameter(mzTabInstrumentSource_)));
    }
    if (mzTabInstrumentAnalyzer_!=null){
      propertyRows.add(new Pair<String,String>(MZTAB_MSANALYZER,getStringFromMzTabParameter(mzTabInstrumentAnalyzer_)));
    }
    if (mzTabInstrumentDetector_!=null){
      propertyRows.add(new Pair<String,String>(MZTAB_DETECTOR,getStringFromMzTabParameter(mzTabInstrumentDetector_)));
    }
    if (alexTargetlist_==true){
      propertyRows.add(new Pair<String,String>(ALEX_TARGETLIST,String.valueOf(alexTargetlist_)));
    }
    if (useMsconvertForWaters_==true) {
      propertyRows.add(new Pair<String,String>(USE_MSCONVERT_FOR_WATERS,String.valueOf(useMsconvertForWaters_)));
    }
    
    String key;
    String value;
    if (faHydroxyEncoding!=null) {
      for (Short oh : faHydroxyEncoding.getHydroxyNumbersInAscendingOrder()) {
        try {
          key = EXCEL_HYDROXY_FA_PREFIX+String.valueOf(oh);
          value = faHydroxyEncoding.getEncodedPrefix(oh);
          propertyRows.add(new Pair<String,String>(key,value));
        //this catch can never happen
        }catch (HydroxylationEncodingException e) {}
      }
    }
    if (lcbHydroxyEncoding!=null) {
      for (Short oh : lcbHydroxyEncoding.getHydroxyNumbersInAscendingOrder()) {
        try {
          key = EXCEL_HYDROXY_LCB_PREFIX+String.valueOf(oh);
          value = lcbHydroxyEncoding.getEncodedPrefix(oh);
          propertyRows.add(new Pair<String,String>(key,value));
        //this catch can never happen
        }catch (HydroxylationEncodingException e) {}
      }
    }
    return propertyRows;
  }
  
  
  /**
   * creates a LipidomicsConstants object from an Excel sheet - this cannot be done in a separate
   * class, since the parameters to be set are private
   * @param sheet Excel sheet to read the parameters from
   * @throws SettingsException thrown when a settings combination is not possible
   * @return settings: [0] LipidomicsConstants object containing the parameters that were read; [1] FA hydroxylation encoding; [2] LCB hydroxylation encoding
   */
  @Deprecated
  //slower Apache POI reader is only used for old xls excel files
  //TODO: remove completely after an adequate transition period *note written: 25.05.2022*
  //(1 year should suffice as xls files are not written since a few years already and newer LDA versions do not support such old file formats anymore anyway) 
  public static Object[] readSettingsFromExcelApachePOI(org.apache.poi.ss.usermodel.Sheet sheet) throws SettingsException{
    Properties properties = new Properties();
    int keyColumn = -1;
    int valueColumn = -1;
    Hashtable<String,Short> faOhEncondings = new Hashtable<String,Short>();
    Hashtable<String,Short> lcbOhEncondings = new Hashtable<String,Short>();
    String ohNumberString;

    for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
      Hashtable<Integer,Object> entries = ExcelUtils.getEntriesOfOneRow(sheet.getRow(rowCount),true);
      // the columns are known - parse the entries
      if (keyColumn>-1 && valueColumn>-1){
        if (!entries.containsKey(keyColumn) || !entries.containsKey(valueColumn) || !(entries.get(keyColumn) instanceof String)) continue;
        String key = (String)entries.get(keyColumn);
        String value = "";
        if (entries.get(valueColumn) instanceof String) value = (String)entries.get(valueColumn);
        else value = String.valueOf((Double)entries.get(valueColumn));
        if (key.startsWith(EXCEL_HYDROXY_FA_PREFIX) || key.startsWith(EXCEL_HYDROXY_LCB_PREFIX)) {
          ohNumberString = null;
          if (key.startsWith(EXCEL_HYDROXY_FA_PREFIX))
            ohNumberString = key.substring(EXCEL_HYDROXY_FA_PREFIX.length());
          else if (key.startsWith(EXCEL_HYDROXY_LCB_PREFIX))
            ohNumberString = key.substring(EXCEL_HYDROXY_LCB_PREFIX.length());
          Short ohNumber = Short.parseShort(ohNumberString);
          if (key.startsWith(EXCEL_HYDROXY_FA_PREFIX))
            faOhEncondings.put(value, ohNumber);
          else if (key.startsWith(EXCEL_HYDROXY_LCB_PREFIX))
            lcbOhEncondings.put(value, ohNumber);
        }else
          properties.put(key, value);
      // looking for the header columns
      } else{
        for (Integer columnId : entries.keySet()){
          if (entries.get(columnId)==null || (!(entries.get(columnId) instanceof String)) || ((String)entries.get(columnId)).length()==0) continue;
          String entry = (String)entries.get(columnId);
          if (entry.equalsIgnoreCase(EXCEL_KEY)) keyColumn = columnId;
          if (entry.equalsIgnoreCase(EXCEL_VALUE)) valueColumn = columnId;
        }
      }
    }
    LipidomicsConstants consts = new LipidomicsConstants(false);
    consts.ldaVersion_ = properties.getProperty(LDA_VERSION, Settings.VERSION);
    String rawFile = null;
    if (properties.containsKey(RAW_FILE)) rawFile = properties.getProperty(RAW_FILE);
    consts.rawFileName_ = rawFile;
    String usedCutoff = null;
    if (properties.containsKey(BASE_PEAK_CUTOFF)) usedCutoff = properties.getProperty(BASE_PEAK_CUTOFF);
    if (usedCutoff!=null){
      try{
        Double.parseDouble(usedCutoff);
        consts.relativeMS1BasePeakCutoff_ = usedCutoff; 
      }catch (NumberFormatException nfx){}
    }
    if (properties.containsKey(MASS_SHIFT))properties.put(MASS_SHIFT, String.valueOf(Double.parseDouble(properties.getProperty(MASS_SHIFT))*1000d));
    consts.setVariables(properties);
    Object[] returnValues = new Object[3];
    returnValues[0] = consts;
    returnValues[1] = new HydroxyEncoding(faOhEncondings);
    returnValues[2] = new HydroxyEncoding(lcbOhEncondings);
    return returnValues;
  }

  
  /**
   * Automatically generated by Eclipse. 
   * Compares all dynamic fields (except those generated by reading the mzTabConfFile) of this class with another Object.
   * The fields generated by reading the mzTabConfFile are excluded, because this method aims at comparing Quantification Result files. The mzTabConfFile is a separate file.
   */
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    LipidomicsConstants other = (LipidomicsConstants) obj;
    return Objects.equals(alexTargetlistUsed_, other.alexTargetlistUsed_)
        && alexTargetlist_ == other.alexTargetlist_
        && Objects.equals(basePeakDefaultCutoff_, other.basePeakDefaultCutoff_)
        && Float.floatToIntBits(broadChromIntensityCutoffNoSmall_) == Float
            .floatToIntBits(other.broadChromIntensityCutoffNoSmall_)
        && Float.floatToIntBits(broadChromIntensityCutoff_) == Float
            .floatToIntBits(other.broadChromIntensityCutoff_)
        && broadChromMeanSmoothRepeats_ == other.broadChromMeanSmoothRepeats_
        && broadChromSmoothRange_ == other.broadChromSmoothRange_
        && broadChromSmoothRepeats_ == other.broadChromSmoothRepeats_
        && Float.floatToIntBits(broadChromSteepnessChangeNoSmall_) == Float
            .floatToIntBits(other.broadChromSteepnessChangeNoSmall_)
        && Float.floatToIntBits(broaderProfileTimeTolerance_) == Float
            .floatToIntBits(other.broaderProfileTimeTolerance_)
        && Double.doubleToLongBits(chainCutoffValue_) == Double
            .doubleToLongBits(other.chainCutoffValue_)
        && checkChainLabelCombination_ == other.checkChainLabelCombination_
        && chromExportShowLegend_ == other.chromExportShowLegend_
        && chromLowestResolution_ == other.chromLowestResolution_
        && chromMultiplicationFactorForInt_ == other.chromMultiplicationFactorForInt_
        && Float.floatToIntBits(chromSmoothRange_) == Float
            .floatToIntBits(other.chromSmoothRange_)
        && chromSmoothRepeats_ == other.chromSmoothRepeats_
        && closePeakTimeTolerance_ == other.closePeakTimeTolerance_
        && Float.floatToIntBits(coarseChromMzTolerance_) == Float
            .floatToIntBits(other.coarseChromMzTolerance_)
        && Objects.equals(currentMSMachine_, other.currentMSMachine_)
        && cvBtoRequired_ == other.cvBtoRequired_
        && cvChmoRequired_ == other.cvChmoRequired_
        && cvClRequired_ == other.cvClRequired_
        && cvDoidRequired_ == other.cvDoidRequired_
        && cvNcbiTaxonRequired_ == other.cvNcbiTaxonRequired_
        && cvSepRequired_ == other.cvSepRequired_
        && Float.floatToIntBits(finalProbeMzCompTolerance_) == Float
            .floatToIntBits(other.finalProbeMzCompTolerance_)
        && Float.floatToIntBits(finalProbeTimeCompTolerance_) == Float
            .floatToIntBits(other.finalProbeTimeCompTolerance_)
        && Objects.equals(intermediateFileFormat_,
            other.intermediateFileFormat_)
        && Float.floatToIntBits(isoInBetweenAreaFactor_) == Float
            .floatToIntBits(other.isoInBetweenAreaFactor_)
        && isoInBetweenMaxTimeDistance_ == other.isoInBetweenMaxTimeDistance_
        && isoNearNormalProbeTime_ == other.isoNearNormalProbeTime_
        && isotopeCorrection_ == other.isotopeCorrection_
        && isotopeInBetweenTime_ == other.isotopeInBetweenTime_
        && Objects.equals(ldaVersion_, other.ldaVersion_)
        && Double.doubleToLongBits(massShift_) == Double
            .doubleToLongBits(other.massShift_)
        && maxFileSizeForChromTranslationAtOnceInMB_ == other.maxFileSizeForChromTranslationAtOnceInMB_
        && Objects.equals(minimumRelativeIntensity_,
            other.minimumRelativeIntensity_)
        && minimumThresholdForHighConfidenceRTMatch_ == other.minimumThresholdForHighConfidenceRTMatch_
        && ms2ChromMultiplicationFactorForInt_ == other.ms2ChromMultiplicationFactorForInt_
        && Float.floatToIntBits(ms2IsobarSCExclusionRatio_) == Float
            .floatToIntBits(other.ms2IsobarSCExclusionRatio_)
        && Float.floatToIntBits(ms2IsobarSCFarExclusionRatio_) == Float
            .floatToIntBits(other.ms2IsobarSCFarExclusionRatio_)
        && Float.floatToIntBits(ms2IsobaricOtherRtDifference_) == Float
            .floatToIntBits(other.ms2IsobaricOtherRtDifference_)
        && ms2MinIntsForNoiseRemoval_ == other.ms2MinIntsForNoiseRemoval_
        && Float.floatToIntBits(ms2MzTolerance_) == Float
            .floatToIntBits(other.ms2MzTolerance_)
        && Float.floatToIntBits(ms2PrecursorTolerance_) == Float
            .floatToIntBits(other.ms2PrecursorTolerance_)
        && ms2_ == other.ms2_
        && Objects.equals(mzTabInstrumentAnalyzer_,
            other.mzTabInstrumentAnalyzer_)
        && Objects.equals(mzTabInstrumentDetector_,
            other.mzTabInstrumentDetector_)
        && Objects.equals(mzTabInstrumentName_, other.mzTabInstrumentName_)
        && Objects.equals(mzTabInstrumentSource_, other.mzTabInstrumentSource_)
        && Objects.equals(mzTabSample_, other.mzTabSample_)
        && Objects.equals(mzUnit_, other.mzUnit_)
        && Float.floatToIntBits(neutronMass_) == Float
            .floatToIntBits(other.neutronMass_)
        && Float.floatToIntBits(noiseCutoffDeviationValue_) == Float
            .floatToIntBits(other.noiseCutoffDeviationValue_)
        && Float.floatToIntBits(overlapDistanceDeviationFactor_) == Float
            .floatToIntBits(other.overlapDistanceDeviationFactor_)
        && Float.floatToIntBits(overlapFullDistanceDivisor_) == Float
            .floatToIntBits(other.overlapFullDistanceDivisor_)
        && Float.floatToIntBits(overlapPeakDistanceDivisor_) == Float
            .floatToIntBits(other.overlapPeakDistanceDivisor_)
        && Float.floatToIntBits(overlapPossibleIntensityThreshold_) == Float
            .floatToIntBits(other.overlapPossibleIntensityThreshold_)
        && Float.floatToIntBits(overlapSureIntensityThreshold_) == Float
            .floatToIntBits(other.overlapSureIntensityThreshold_)
        && peakDiscardingAreaFactor_ == other.peakDiscardingAreaFactor_
        && Float.floatToIntBits(profileGeneralIntCutoff_) == Float
            .floatToIntBits(other.profileGeneralIntCutoff_)
        && Float.floatToIntBits(profileIntThreshold_) == Float
            .floatToIntBits(other.profileIntThreshold_)
        && Float.floatToIntBits(profileIntensityCutoff1_) == Float
            .floatToIntBits(other.profileIntensityCutoff1_)
        && Float.floatToIntBits(profileIntensityCutoff2_) == Float
            .floatToIntBits(other.profileIntensityCutoff2_)
        && Float.floatToIntBits(profileMaxRange_) == Float
            .floatToIntBits(other.profileMaxRange_)
        && profileMeanSmoothRepeats_ == other.profileMeanSmoothRepeats_
        && Float.floatToIntBits(profileMzMinRange_) == Float
            .floatToIntBits(other.profileMzMinRange_)
        && Float.floatToIntBits(profileMzRange_) == Float
            .floatToIntBits(other.profileMzRange_)
        && Float.floatToIntBits(profilePeakAcceptanceRange_) == Float
            .floatToIntBits(other.profilePeakAcceptanceRange_)
        && Float.floatToIntBits(profileSmoothRange_) == Float
            .floatToIntBits(other.profileSmoothRange_)
        && profileSmoothRepeats_ == other.profileSmoothRepeats_
        && Float.floatToIntBits(profileSmoothingCorrection_) == Float
            .floatToIntBits(other.profileSmoothingCorrection_)
        && Float.floatToIntBits(profileSteepnessChange1_) == Float
            .floatToIntBits(other.profileSteepnessChange1_)
        && Float.floatToIntBits(profileSteepnessChange2_) == Float
            .floatToIntBits(other.profileSteepnessChange2_)
        && Float.floatToIntBits(profileTimeTolerance_) == Float
            .floatToIntBits(other.profileTimeTolerance_)
        && Objects.equals(rawFileName_, other.rawFileName_)
        && Float.floatToIntBits(relativeAreaCutoff_) == Float
            .floatToIntBits(other.relativeAreaCutoff_)
        && Float.floatToIntBits(relativeFarAreaCutoff_) == Float
            .floatToIntBits(other.relativeFarAreaCutoff_)
        && relativeFarAreaTimeSpace_ == other.relativeFarAreaTimeSpace_
        && Float.floatToIntBits(relativeIsoInBetweenCutoff_) == Float
            .floatToIntBits(other.relativeIsoInBetweenCutoff_)
        && Objects.equals(relativeMS1BasePeakCutoff_,
            other.relativeMS1BasePeakCutoff_)
        && removeFromOtherIsotopes_ == other.removeFromOtherIsotopes_
        && respectIsotopicDistribution_ == other.respectIsotopicDistribution_
        && scanStep_ == other.scanStep_
        && shotgunIntensityRemoval_ == other.shotgunIntensityRemoval_
        && shotgunProcessing_ == other.shotgunProcessing_
        && Float.floatToIntBits(shotgunRelIntCutoff_) == Float
            .floatToIntBits(other.shotgunRelIntCutoff_)
        && shotgun_ == other.shotgun_
        && Float.floatToIntBits(smallChromIntensityCutoff_) == Float
            .floatToIntBits(other.smallChromIntensityCutoff_)
        && smallChromMeanSmoothRepeats_ == other.smallChromMeanSmoothRepeats_
        && Float.floatToIntBits(smallChromMzRange_) == Float
            .floatToIntBits(other.smallChromMzRange_)
        && Float.floatToIntBits(smallChromSmoothRange_) == Float
            .floatToIntBits(other.smallChromSmoothRange_)
        && smallChromSmoothRepeats_ == other.smallChromSmoothRepeats_
        && sparseData_ == other.sparseData_
        && Objects.equals(threeDViewerDefaultMZResolution_,
            other.threeDViewerDefaultMZResolution_)
        && Objects.equals(threeDViewerDefaultTimeResolution_,
            other.threeDViewerDefaultTimeResolution_)
        && Objects.equals(threeDViewerMs2DefaultMZResolution_,
            other.threeDViewerMs2DefaultMZResolution_)
        && Objects.equals(threeDViewerMs2DefaultTimeResolution_,
            other.threeDViewerMs2DefaultTimeResolution_)
        && Float.floatToIntBits(twinInBetweenCutoff_) == Float
            .floatToIntBits(other.twinInBetweenCutoff_)
        && Float.floatToIntBits(unionInBetweenCutoff_) == Float
            .floatToIntBits(other.unionInBetweenCutoff_)
        && use3D_ == other.use3D_
        && useMostOverlappingIsotopeOnly_ == other.useMostOverlappingIsotopeOnly_
        && useMsconvertForWaters_ == other.useMsconvertForWaters_
        && useNoiseCutoff_ == other.useNoiseCutoff_;
  }

  /**
   * returns the String representation of the Param VO
   * @param param the Param VO
   * @return the String representation of the VO
   */
  private String getStringFromMzTabParameter(Parameter param){
    String returnString = param.getCvLabel()+",";
    if (param.getCvAccession()!=null && param.getCvAccession().length()>0)
      returnString += param.getCvAccession()+",";
    if (param.getName()!=null && param.getName().length()>0)
      returnString += param.getName()+",";
    if (param.getValue()!=null && param.getValue().length()>0)
      returnString += param.getValue()+",";
    return returnString;
  }

  public String getRelativeMS1BasePeakCutoff()
  {
    return relativeMS1BasePeakCutoff_;
  }

  public void setRelativeMS1BasePeakCutoff(String relativeMS1BasePeakCutoff)
  {
    this.relativeMS1BasePeakCutoff_ = relativeMS1BasePeakCutoff;
  }
  
  public void setLDAVersion(String ldaVersion) {
    this.ldaVersion_ = ldaVersion;
  }
  
  public String getLDAVersion() {
	  return ldaVersion_;
  }
  
  public String getRawFileName(){
    return rawFileName_;
  }
  
  public void setRawFileName(String rawFileName){
    rawFileName_ = rawFileName;
  }

  public String getMSMachine()
  {
    return currentMSMachine_;
  }
  
  public double getShift()
  {
    return this.massShift_;
  }
  
    /** was this file quantified by using a Alex123 target list*/
  public boolean isAlexTargetlist()
  {
    return alexTargetlist_;
  }

  /** set whether this file was quantified by using a Alex123 target list*/
  public void setAlexTargetlist(boolean alexTargetlist)
  {
    this.alexTargetlist_ = alexTargetlist;
  }
  
  /** typically, for Waters files Mass++ is used; however, when msconvert shall be used, this returns true*/
  public static boolean useMsconvertForWaters()
  {
    getInstance();
    return instance_.useMsconvertForWaters_;
  }

  /** lookup for classes whether MSn fragments were defined*/
  public Hashtable<String,Boolean> getAlexTargetlistUsed()
  {
    return alexTargetlistUsed_;
  }

  /** set lookup for classes whether MSn fragments were defined*/
  public void setAlexTargetlistUsed(Hashtable<String,Boolean> alexTargetlistUsed)
  {
    this.alexTargetlistUsed_ = alexTargetlistUsed;
  }
  
  /**
   *  returns whether for this search the isotopic distribution has to match for a valid hit
   */
  public boolean getRespectIsotopicDistribution(){
    return respectIsotopicDistribution_;
  }
  
  /**
   * returns the ambigious position separator depending on the LDA version
   */
  public String getChainCombiSeparatorAmbigPosDependingOnVersion()
  {
	  if (shouldOldEncodingBeUsed())
		  return CHAIN_COMBI_SEPARATOR_AMBIG_POS_OLD;
	  
	  // return null for new encoding since it is not a separator anymore (sn nomenclature)
	  return null;	  
  }
  
  /**
   * returns the modification separator depending on the LDA version
   */
  public String getModSeparatorDependingOnVersion()
  {
	  if (!shouldOldEncodingBeUsed())
		  return CHAIN_MOD_SEPARATOR;
	  
	  // old separator (at least starting one)
	  return "[";	  
  }
  
  /**
   * returns true if old encoding should be used, false otherwise
   */
  public boolean shouldOldEncodingBeUsed()
  {
	  String[] version = this.ldaVersion_.split("_")[0].split("\\.");

	  // for versions >= 2.9 use new encoding
	  if (Integer.parseInt(version[0]) > 2)
		  return false;
	  else if (Integer.parseInt(version[0]) >= 2) {
		  if (Integer.parseInt(version[1]) >= 9)
				  return false;
	  }

	  // for all other versions use old encoding
	  return true;
  }
  
  /**
   * 
   * @param tolerance the tolerance value in its original unit
   * @param mz the m/z value to refer to
   * @return the tolerance value in Da
   */
  private static float getCorrectMzTolerance(float tolerance, String mzUnit, float mz){
    if (mzUnit.equalsIgnoreCase(LipidomicsConstants.MZUNIT_DA))
      return tolerance;
    else if (mzUnit.equalsIgnoreCase(LipidomicsConstants.MZUNIT_PPM))
      return (tolerance*mz)/1000000f;
    else
      return tolerance;
  }
  
}
