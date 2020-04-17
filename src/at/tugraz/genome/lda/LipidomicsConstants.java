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
import java.util.Properties;
import java.util.Vector;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;

import de.isas.mztab2.model.Contact;
import de.isas.mztab2.model.Instrument;
import de.isas.mztab2.model.Parameter;
import de.isas.mztab2.model.Publication;
import de.isas.mztab2.model.PublicationItem;
import de.isas.mztab2.model.PublicationItem.TypeEnum;
import de.isas.mztab2.model.Sample;
import de.isas.mztab2.model.SampleProcessing;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsConstants
{
  public final static String EXCEL_MS_OH = "OH";
  public final static int EXCEL_NO_OH_INFO = -1;
  
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
  
  private final static String EXCEL_HYDROXY_FA_PREFIX = "faOHEncoding_";
  private final static String EXCEL_HYDROXY_LCB_PREFIX = "lcbOHEncoding_";
  
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
  public final static String CHAIN_COMBI_SEPARATOR_AMBIG_POS = ";";
  
  /** the String for separating the OH numbers*/
  public final static String ALEX_OH_SEPARATOR = ";";
  /** prefix for alkylated chains*/
  public final static String ALEX_ALKYL_PREFIX = "O_";
  /** prefix for alkenylated chains*/
  public final static String ALEX_ALKENYL_PREFIX = "P_";
  /** the String for separating the chains*/
  public final static String ALEX_CHAIN_SEPARATOR = "-";

  
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
  private float twinPeakMzTolerance_;
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
  
  
  /** this is the conventional chromatography mode*/
  public final static short SHOTGUN_FALSE = 0;
  /** this is the shotgun mode*/
  public final static short SHOTGUN_TRUE = 1;
  /** this is the PRM mode*/
  public final static short SHOTGUN_PRM = 2;
  
  
  private final static String LDA_VERSION = "LDA-version";
  private final static String RAW_FILE = "rawFile";
  private final static String BASE_PEAK_CUTOFF = "basePeakCutoff";
  private final static String BASE_PEAK_CUTOFF_DEFAULT = "0.1";
  private final static String MASS_SHIFT = "massShift";
  private final static String MASS_SHIFT_DEFAULT = "0";
  private final static String MACHINE_NAME = "machineName";
  private final static String MACHINE_NAME_DEFAULT = null;
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
  private final static String TWIN_PEAK_MZ_TOL = "twinPeakMzTolerance";
  private final static String TWIN_PEAK_MZ_TOL_DEFAULT = "0.002";
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

  
  private final static String MZTAB_INSTRUMENT = "mzTabInstrumentName";
  private final static String MZTAB_IONSOURCE = "mzTabInstrumentIonsource";
  private final static String MZTAB_MSANALYZER = "mzTabInstrumentAnalyzer";
  private final static String MZTAB_DETECTOR = "mzTabInstrumentDetector";
  
  private final static String CHROM_EXPORT_LEGEND_IN_CHROM = "chromexportPaintLegendInChrom";
  
  private final static String EXCEL_KEY = "Key";
  private final static int EXCEL_KEY_COLUMN = 0;
  private final static String EXCEL_VALUE = "Value";
  private final static int EXCEL_VALUE_COLUMN = 1;
  
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
  
  public static LipidomicsConstants getInstance() {
    if (instance_ == null) {
      instance_ = new LipidomicsConstants();
    }
    return instance_;
  }
  
  private LipidomicsConstants(boolean dummy){
    cvSepRequired_ = false;
    cvNcbiTaxonRequired_ = false;
    cvClRequired_ = false;
    cvBtoRequired_ = false;  
    cvDoidRequired_ = false;
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
  
  private void setVariables(Properties properties) throws SettingsException{
    
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
    twinPeakMzTolerance_ = Float.parseFloat(properties.getProperty(TWIN_PEAK_MZ_TOL,TWIN_PEAK_MZ_TOL_DEFAULT));
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
   * @return the file pieces for the file translation to the chrom file
   */
  public static int getmMaxFileSizeForChromTranslationAtOnceInMB()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.maxFileSizeForChromTranslationAtOnceInMB_;
  }
  
  /**
   * @return the m/z tolerance for the first coarse chromatogram
   */
  public static float getCoarseChromMzTolerance(float mz)
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return getCorrectMzTolerance(instance_.coarseChromMzTolerance_,instance_.mzUnit_,mz);
  }
  
  public static double getMassShift()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.massShift_;
  }
  
  public static float getChromSmoothRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.chromSmoothRange_;
  }

  public static int getChromSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.chromSmoothRepeats_;
  }

  public static boolean isotopicCorrection()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.isotopeCorrection_;
  }
  
  public static boolean removeIfOtherIsotopePresent()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.removeFromOtherIsotopes_;
  }

  public static boolean removeIfDistriDoesNotFit()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.respectIsotopicDistribution_;
  }
  
  public static boolean useNoiseCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.useNoiseCutoff_;
  }
  
  /**
   * 
   * @return check if number of labels in species name corresponds with number of labels in chains
   */
  public static boolean checkChainLabelCombination()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.checkChainLabelCombination_;
  }
  
  public static float getNoiseCutoffDeviationValue()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.noiseCutoffDeviationValue_;
  }
  
  public static Float getMinimumRelativeIntensity()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.minimumRelativeIntensity_;
  }
  

  public static int getScanStep()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.scanStep_;
  }
  
  /**
   * @return if 3D-method should be used for automated quantification
   */
  public static boolean use3D()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.use3D_;
  }
  
  /**
   * @return the m/z tolerance for the first coarse chromatogram
   */
  public static float getProfileMzRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileMzRange_;
  }
  
  /**
   * @return the time tolerance of the profile
   */
  public static float getProfileTimeTolerance_()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileTimeTolerance_;
  }
  
  /**
   * @return the intensity threshold for the extraction of a new the profile
   */
  public static float getProfileIntThreshold_()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileIntThreshold_;
  }
  
  /**
   * @return the range for the broader profile extraction
   */
  public static float getBroaderProfileTimeTolerance_()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broaderProfileTimeTolerance_;
  }
  
  /**
   * @return the m/z smoothing range for the profile
   */
  public static float getProfileSmoothRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileSmoothRange_;
  }
  
  /**
   * @return the amount of smoothing repeats for the profile
   */
  public static int getProfileSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileSmoothRepeats_;
  }
  /**
   * @return the amount of smoothing repeats for the profile
   */
  public static int getProfileMeanSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileMeanSmoothRepeats_;
  }  
  /**
   * @return the amount of smoothing repeats for the profile
   */
  public static float getProfileMzMinRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileMzMinRange_;
  }
  /**
   * @return steepness change1 threshold for profile
   */
  public static float getProfileSteepnessChange1()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileSteepnessChange1_;
  }
  
  /**
   * @return steepness change2 threshold for profile
   */
  public static float getProfileSteepnessChange2()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileSteepnessChange2_;
  }  
  
  /**
   * @return steepness intensity cutoff 1 threshold for profile
   */
  public static float getProfileIntensityCutoff1()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileIntensityCutoff1_;
  }  

  /**
   * @return steepness intensity cutoff 2 threshold for profile
   */
  public static float getProfileIntensityCutoff2()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileIntensityCutoff2_;
  }  

  /**
   * @return steepness general intensity cutoff threshold for profile
   */
  public static float getProfileGeneralIntCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileGeneralIntCutoff_;
  }
  
  /**
   * @return steepness general intensity cutoff threshold for profile
   */
  public static float getProfilePeakAcceptanceRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profilePeakAcceptanceRange_;
  }
  
  /**
   * @return m/z correction for profile smoothing
   */
  public static float getProfileSmoothingCorrection()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileSmoothingCorrection_;
  }
  
  
  /**
   * @return the maximum size of one profile
   */
  public static float getProfileMaxRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.profileMaxRange_;
  }
  
  /**
   * @return m/z range for the small chromatogram
   */
  public static float getSmallChromMzRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.smallChromMzRange_;
  } 
  
  /**
   * @return smoothing repeats for the small chromatogram
   */
  public static int getSmallChromSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.smallChromSmoothRepeats_;
  }
  
  /**
   * @return smoothing repeats for the small chromatogram
   */
  public static int getSmallChromMeanSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.smallChromMeanSmoothRepeats_;
  }
  
  /**
   * @return smoothing repeats for the small chromatogram
   */
  public static float getSmallChromSmoothRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.smallChromSmoothRange_;
  } 
  
  /**
   * @return intensity-cutoff for the small chromatogram
   */
  public static float getSmallChromIntensityCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.smallChromIntensityCutoff_;
  }
  
  /**
   * @return smoothing repeats for the broad chromatogram
   */
  public static int getBroadChromSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broadChromSmoothRepeats_;
  }
  
  /**
   * @return smoothing repeats for the broad chromatogram
   */
  public static int getBroadChromMeanSmoothRepeats()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broadChromMeanSmoothRepeats_;
  }
  
  /**
   * @return smoothing range for the broad chromatogram
   */
  public static float getBroadChromSmoothRange()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broadChromSmoothRange_;
  } 
  
  /**
   * @return intensity-cutoff for the broad chromatogram
   */
  public static float getBroadChromIntensityCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broadChromIntensityCutoff_;
  }  
  
  /**
   * @return steepness-change threshold for the broad chromatogram if one of the 
   * small borders is not valid
   */
  public static float getBroadChromSteepnessChangeNoSmall()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broadChromSteepnessChangeNoSmall_;
  }
  
  /**
   * @return intensity-cutoff for the broad chromatogram if one of the 
   * small borders is not valid
   */
  public static float getBroadIntensityCutoffNoSmall()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.broadChromIntensityCutoffNoSmall_;
  }
  
  /**
   * @return the tolerance or a comparison of the final probe in time direction
   */
  public static float getFinalProbeTimeCompTolerance()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.finalProbeTimeCompTolerance_;
  }
  
  /**
   * @return the tolerance for a comparison of the final probe in m/z direction
   */
  public static float getFinalProbeMzCompTolerance()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.finalProbeMzCompTolerance_;
  }

  
  /**
   * @return the distance tolerance factor for the detection of an overlap;
   * the comparison is done between the original probe and the newly calculated one
   */
  public static float getOverlapDistanceDeviationFactor()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.overlapDistanceDeviationFactor_;
  }
  
  /**
   * @return a multiplication factor for the comparison of intensities
   * determines if an overlap is possible
   */
  public static float getOverlapPossibleIntensityThreshold()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.overlapPossibleIntensityThreshold_;
  }

  /**
   * @return a multiplication factor for the comparison of intensities
   * determines if an overlap is sure
   */
  public static float getOverlapSureIntensityThreshold()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.overlapSureIntensityThreshold_;
  }
  /**
   * @return a divisor to check if the peak of an isotopic peak is near
   * the other peak to assume an overlap
   * this distance is between peak and one border
   */
  public static float getOverlapPeakDistanceDivisor()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.overlapPeakDistanceDivisor_;
  }
  /**
   * @return a divisor to check if the peak of an isotopic peak is near
   * the other peak to assume an overlap
   * this distance is between the borders
   */
  public static float getOverlapFullDistanceDivisor()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.overlapFullDistanceDivisor_;
  }
  
  /**
   * @return peaks that are relatively smaller to the highest found peak than this factor are discarded
   */
  public static int getPeakDiscardingAreaFactor()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.peakDiscardingAreaFactor_;
  }
  /**
   * @return allowed time distance in seconds to check if there is an isotope between the two peaks
   */
  public static int getIsotopeInBetweenTime()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.isotopeInBetweenTime_;
  }
  /**
   * @return if an isotope is in between the smaller one of the two peaks is regarded to belong to be a fragment
   * of the other isotope if the area of the smaller one is isoInBetweenAreaFactor times smaller than
   * the bigger one
   */
  public static float getIsoInBetweenAreaFactor()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.isoInBetweenAreaFactor_;
  }
  
  public static int getIsoInBetweenMaxTimeDistance()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.isoInBetweenMaxTimeDistance_;
  }
  
  /**
   * @return this is the time distance in seconds where an isotopic peak is still regarded as interesting
   */
  public static int getIsoNearNormalProbeTime()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.isoNearNormalProbeTime_;
  }
  /**
   * @return after the isotopic peaks are thrown out a higher threshold is applied to discard small peaks 
   * (relative to the highest found one)
   */
  public static float getRelativeAreaCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.relativeAreaCutoff_;
  }
  /**
   * @return intensities that are close to the highest one can be members of a twin peak, hence they are
   * not discarded so quickly
   * this is a relative threshold to discard peaks that are farer away
   */
  public static float getRelativeFarAreaCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.relativeFarAreaCutoff_;
  }
  /**
   * @return intensities that are close to the highest one can be members of a twin peak, hence they are
   * not discarded so quickly
   * this is the time distance in seconds to define if a peak is near the highest one
   */
  public static int getRelativeFarAreaTimeSpace()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.relativeFarAreaTimeSpace_;
  }
  public static float getTwinInBetweenCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.twinInBetweenCutoff_;
  }
  public static float getUnionInBetweenCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.unionInBetweenCutoff_;
  }
  /**
   * @return if still there is an isotopic peak in between the smaller peak is normally discarded
   * if the ratio between smallerPeak/higherPeak the smaller one is not discarded since
   * it cannot be clearly identified which one of the peaks is the correct one
   */
  public static float getRelativeIsoInBetweenCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.relativeIsoInBetweenCutoff_;
  }

  /**
   * @return tolerance in m/z direction to allow the union of twin peaks
   */
  public static float getTwinPeakMzTolerance()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.twinPeakMzTolerance_;
  }
  /**
   * @return this is the time distance in seconds defines if in a last step close peaks should be united
   */
  public static int getClosePeakTimeTolerance()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.closePeakTimeTolerance_;
  }
  
  public static String getBasePeakDefaultCutoff()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.basePeakDefaultCutoff_;
  }
  public static int getChromMultiplicationFactorForInt(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.chromMultiplicationFactorForInt_;
  }

  public static int getChromLowestResolution(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.chromLowestResolution_;
  }

  public static String getThreeDViewerDefaultTimeResolution(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.threeDViewerDefaultTimeResolution_;
  }
  
  public static String getThreeDViewerDefaultMZResolution(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.threeDViewerDefaultMZResolution_;
  }
  
  /**
   * @return if MS2 should be used for analyisis
   */
  public static boolean isMS2()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2_;
  }
  
  public static float getMs2PrecursorTolerance(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2PrecursorTolerance_;
  }
  
  public static int getMs2ChromMultiplicationFactorForInt(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2ChromMultiplicationFactorForInt_;
  }
  
  public static String getThreeDViewerMs2DefaultTimeResolution_(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.threeDViewerMs2DefaultTimeResolution_;
  }

  public static String getThreeDViewerMs2DefaultMZResolution(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.threeDViewerMs2DefaultMZResolution_;
  }
  
  public static boolean isChromExportShowLegend()
  {
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.chromExportShowLegend_;
  }
  
  /**
   * 
   * @return the +/- m/z tolerance for peak detection
   */
  public static float getMs2MzTolerance(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2MzTolerance_;
  }
  
  /**
   * 
   * @return minimum number of detected signals to start a noise removal
   */
  public static int getMs2MinIntsForNoiseRemoval(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2MinIntsForNoiseRemoval_;
  }
  
  /**
   * 
   * @return relative intensity cutoff for exclusion of isobars regarding spectrum coverage
   */
  public static float getMs2IsobarSCExclusionRatio(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2IsobarSCExclusionRatio_;
  }

  
  /**
   * 
   * @return relative intensity cutoff for exclusion of isobars regarding spectrum coverage, for areas that are farer away from a unique peak identification
   */
  public static float getMs2IsobarSCFarExclusionRatio(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2IsobarSCFarExclusionRatio_;
  }

  /**
   * 
   * @return retention time in minutes that define a peak to be farer away, so that the ms2IsobarSCFarExclusionRatio can be used
   */
  public static float getMs2IsobaricOtherRtDifference(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.ms2IsobaricOtherRtDifference_;
  }
  
  /**
   * 
   * @return 0 when LC, 1 when shotgun, 2 when PRM
   */
  public static short isShotgun(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.shotgun_;
  }
  
  /**
   * 
   * @return the unit type of the m/z value
   */
  public static String getMzUnit(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.mzUnit_;    
  }
  
  /**
   * 
   * @return the shotgun processing type according to the definitions in LipidomicsAnalyzer
   */
  public static int getShogunProcessing(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.shotgunProcessing_;    
  }
  
  /**
   * 
   * @return should shotgun intensities be removed below a certain threshold
   */
  public static boolean isShotgunIntensityRemoval(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.shotgunIntensityRemoval_;    
  }

  /**
   * 
   * @return the relative intensity cutoff for the removal
   */
  public static float getShotgunRelIntCutoff(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.shotgunRelIntCutoff_;    
  }

  /**
   * 
   * @return relative cutoff threshold for fatty acid chain detection - it is in relation to the most intense chain
   */
  public static double getChainCutoffValue(){
    if (instance_ == null) LipidomicsConstants.getInstance();
    return instance_.chainCutoffValue_;
  }
  
  public static String getCurrentMSMachine(){
    LipidomicsConstants.getInstance();
    return instance_.currentMSMachine_;
  }
  
  public static float getNeutronMass(){
    LipidomicsConstants.getInstance();
    return instance_.neutronMass_;
  }
  
  public static boolean useMostOverlappingIsotopeOnly(){
    LipidomicsConstants.getInstance();    
    return instance_.useMostOverlappingIsotopeOnly_;
  }
  
  /**
   * 
   * @return contains the acquired data only sparse MS1 data points
   */
  public static boolean isSparseData(){
    LipidomicsConstants.getInstance();
    return instance_.sparseData_;
  }
  
  public static Instrument getMzTabInstrument(){
    LipidomicsConstants.getInstance();
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
    LipidomicsConstants.getInstance();
    return instance_.mzTabContacts_;
  }

  /**
   * 
   * @return for mzTab-export: list of sample processing steps defined in the MZTAB_PROPERTIES_FILE
   */
  public static List<SampleProcessing> getMzTabSampleprocessings(){
    LipidomicsConstants.getInstance();
    return instance_.mzTabSampleProcessings_;
  }

  /**
   * 
   * @return for mzTab-export: list of publications this data refers to
   */
  public static List<Publication> getMzTabPublications(){
    LipidomicsConstants.getInstance();
    return instance_.mzTabPubs_;
  }
  
  /**
   * 
   * @return for mzTab-export: list of fragmentation methods used for the data
   */
  public static List<Parameter> getFragmentationMethods(){
    LipidomicsConstants.getInstance();
    return instance_.fragmethods_;
  }

  /**
   * 
   * @return for mzTab-export: the sample this data originates of
   */
  public static Sample getMzTabSample(){
    LipidomicsConstants.getInstance();
    return instance_.mzTabSample_;    
  }
  
  /**
   * 
   * @return true when this is a shotgun instance
   */
  public short getShotgun(){
    LipidomicsConstants.getInstance();
    return this.shotgun_;
  }
  
  /**
   * checks the adduct lookup and returns the corresponding mzTab notation when present, null otherwise
   * @param ldaAdduct the LDA adduct notation
   * @return the mzTab notation
   */
  public static String getMzTabAdduct(String ldaAdduct){
    LipidomicsConstants.getInstance();
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
    LipidomicsConstants.getInstance();
    return instance_.cvSepRequired_;    
  }

  /**
   * 
   * @return contains CV for NCBITaxon ontology
   */
  public static boolean isCvNcbiTaxonRequired(){
    LipidomicsConstants.getInstance();
    return instance_.cvNcbiTaxonRequired_;
  }

  /**
   * 
   * @return contains CV for CL ontology
   */
  public static boolean isClRequired(){
    LipidomicsConstants.getInstance();
    return instance_.cvClRequired_;
  }
  
  /**
   * 
   * @return contains CV for BTO ontology
   */
  public static boolean isCvBtoRequired(){
    LipidomicsConstants.getInstance();
    return instance_.cvBtoRequired_;
  }

  /**
   * 
   * @return contains CV for DOID ontology
   */
  public static boolean isCvDoidRequired(){
    LipidomicsConstants.getInstance();
    return instance_.cvDoidRequired_;
  }
  
  public static void switchToOtherConfFile(File newConfFile){
    instance_.readConstantsFile(newConfFile.getAbsolutePath());
  }
    
  public static void switchToOtherDefaultConfFile(File newConfFile) throws IOException{
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
  
  /**
   * writes the current settings into an Excel sheet - this cannot be done in a separate
   * class, since the parameters to write are private
   * @param sheet the Excel sheet to be written
   * @param headerStyle style for the header column
   * @param faHydroxyEncoding the character encoding of the number of hydroxylation sites for the
   * @param lcbHydroxyEncoding the character encoding of the number of hydroxylation sites for the LCB
   */
  public void writeSettingsToExcel(Sheet sheet, CellStyle headerStyle, HydroxyEncoding faHydroxyEncoding,
      HydroxyEncoding lcbHydroxyEncoding){
    int rowCount = 0;
    Row row = sheet.createRow(rowCount);
    rowCount++;
    Cell cell = row.createCell(EXCEL_KEY_COLUMN);
    cell.setCellValue(EXCEL_KEY);
    cell.setCellStyle(headerStyle);
    cell = row.createCell(EXCEL_VALUE_COLUMN);
    cell.setCellValue(EXCEL_VALUE);
    cell.setCellStyle(headerStyle);
    int longestKey = 0;
    int longestValue = 0;
    
    rowCount = createPropertyRow(sheet,rowCount,LDA_VERSION,ldaVersion_);
    if (LDA_VERSION.length()>longestKey) longestKey = LDA_VERSION.length();
    if (ldaVersion_.length()>longestValue) longestValue = ldaVersion_.length();
    if (rawFileName_!=null){
      rowCount = createPropertyRow(sheet,rowCount,RAW_FILE,rawFileName_);
      if (RAW_FILE.length()>longestKey) longestKey = RAW_FILE.length();
      if (rawFileName_.length()>longestValue) longestValue = rawFileName_.length();
    }
    rowCount = createPropertyRow(sheet,rowCount,MACHINE_NAME,currentMSMachine_);
    if (MACHINE_NAME.length()>longestKey) longestKey = MACHINE_NAME.length();
    if (currentMSMachine_.length()>longestValue) longestValue = currentMSMachine_.length();
    rowCount = createPropertyRow(sheet,rowCount,NEUTRON_MASS,String.valueOf(neutronMass_));
    if (NEUTRON_MASS.length()>longestKey) longestKey = NEUTRON_MASS.length();
    if (String.valueOf(neutronMass_).length()>longestValue) longestValue = String.valueOf(neutronMass_).length();
    rowCount = createPropertyRow(sheet,rowCount,COARSE_CHROM_MZ_TOL,String.valueOf(coarseChromMzTolerance_));
    if (COARSE_CHROM_MZ_TOL.length()>longestKey) longestKey = COARSE_CHROM_MZ_TOL.length();
    if (String.valueOf(coarseChromMzTolerance_).length()>longestValue) longestValue = String.valueOf(coarseChromMzTolerance_).length();
    rowCount = createPropertyRow(sheet,rowCount,MS2,String.valueOf(ms2_));
    if (MS2.length()>longestKey) longestKey = MS2.length();
    if (String.valueOf(ms2_).length()>longestValue) longestValue = String.valueOf(ms2_).length();
    if (relativeMS1BasePeakCutoff_!=null){
      rowCount = createPropertyRow(sheet,rowCount,BASE_PEAK_CUTOFF,String.valueOf(relativeMS1BasePeakCutoff_));
      if (BASE_PEAK_CUTOFF.length()>longestKey) longestKey = BASE_PEAK_CUTOFF.length();
    }
    if (String.valueOf(basePeakDefaultCutoff_).length()>longestValue) longestValue = String.valueOf(basePeakDefaultCutoff_).length();
    rowCount = createPropertyRow(sheet,rowCount,MASS_SHIFT,String.valueOf(massShift_));
    if (MASS_SHIFT.length()>longestKey) longestKey = MASS_SHIFT.length();
    if (String.valueOf(massShift_).length()>longestValue) longestValue = String.valueOf(massShift_).length();
    rowCount = createPropertyRow(sheet,rowCount,VIEWER_TIME_RESOLUTION,String.valueOf(threeDViewerDefaultTimeResolution_));
    if (VIEWER_TIME_RESOLUTION.length()>longestKey) longestKey = VIEWER_TIME_RESOLUTION.length();
    if (String.valueOf(threeDViewerDefaultTimeResolution_).length()>longestValue) longestValue = String.valueOf(threeDViewerDefaultTimeResolution_).length();
    rowCount = createPropertyRow(sheet,rowCount,VIEWER_MZ_RESOLUTION,String.valueOf(threeDViewerDefaultMZResolution_));
    if (VIEWER_MZ_RESOLUTION.length()>longestKey) longestKey = VIEWER_MZ_RESOLUTION.length();
    if (String.valueOf(threeDViewerDefaultMZResolution_).length()>longestValue) longestValue = String.valueOf(threeDViewerDefaultMZResolution_).length();    
    rowCount = createPropertyRow(sheet,rowCount,MS2_PRECURSOR_TOL,String.valueOf(ms2PrecursorTolerance_));
    if (MS2_PRECURSOR_TOL.length()>longestKey) longestKey = MS2_PRECURSOR_TOL.length();
    if (String.valueOf(ms2PrecursorTolerance_).length()>longestValue) longestValue = String.valueOf(ms2PrecursorTolerance_).length();    
    rowCount = createPropertyRow(sheet,rowCount,MS2_MZ_TOL,String.valueOf(ms2MzTolerance_));
    if (MS2_MZ_TOL.length()>longestKey) longestKey = MS2_MZ_TOL.length();
    if (String.valueOf(ms2MzTolerance_).length()>longestValue) longestValue = String.valueOf(ms2MzTolerance_).length();
    rowCount = createPropertyRow(sheet,rowCount,MS2_MIN_NOISE_REMOVAL,String.valueOf(ms2MinIntsForNoiseRemoval_));
    if (MS2_MIN_NOISE_REMOVAL.length()>longestKey) longestKey = MS2_MIN_NOISE_REMOVAL.length();
    if (String.valueOf(ms2MinIntsForNoiseRemoval_).length()>longestValue) longestValue = String.valueOf(ms2MinIntsForNoiseRemoval_).length();   
    rowCount = createPropertyRow(sheet,rowCount,MS2_ISOBAR_RATIO,String.valueOf(ms2IsobarSCExclusionRatio_));
    if (MS2_ISOBAR_RATIO.length()>longestKey) longestKey = MS2_ISOBAR_RATIO.length();
    if (String.valueOf(ms2IsobarSCExclusionRatio_).length()>longestValue) longestValue = String.valueOf(ms2IsobarSCExclusionRatio_).length();
    rowCount = createPropertyRow(sheet,rowCount,MS2_ISOBAR_FAR_RATIO,String.valueOf(ms2IsobarSCFarExclusionRatio_));
    if (MS2_ISOBAR_FAR_RATIO.length()>longestKey) longestKey = MS2_ISOBAR_FAR_RATIO.length();
    if (String.valueOf(ms2IsobarSCFarExclusionRatio_).length()>longestValue) longestValue = String.valueOf(ms2IsobarSCFarExclusionRatio_).length();
    rowCount = createPropertyRow(sheet,rowCount,MS2_ISOBAR_FAR_RT,String.valueOf(ms2IsobaricOtherRtDifference_));
    if (MS2_ISOBAR_FAR_RT.length()>longestKey) longestKey = MS2_ISOBAR_FAR_RT.length();
    if (String.valueOf(ms2IsobaricOtherRtDifference_).length()>longestValue) longestValue = String.valueOf(ms2IsobaricOtherRtDifference_).length();   
    rowCount = createPropertyRow(sheet,rowCount,CHAIN_CUTOFF,String.valueOf(chainCutoffValue_));
    if (CHAIN_CUTOFF.length()>longestKey) longestKey = CHAIN_CUTOFF.length();
    if (String.valueOf(chainCutoffValue_).length()>longestValue) longestValue = String.valueOf(chainCutoffValue_).length();        
    rowCount = createPropertyRow(sheet,rowCount,MS2_CHROM_MULT_FOR_INT,String.valueOf(ms2ChromMultiplicationFactorForInt_));
    if (MS2_CHROM_MULT_FOR_INT.length()>longestKey) longestKey = MS2_CHROM_MULT_FOR_INT.length();
    if (String.valueOf(ms2ChromMultiplicationFactorForInt_).length()>longestValue) longestValue = String.valueOf(ms2ChromMultiplicationFactorForInt_).length();
    rowCount = createPropertyRow(sheet,rowCount,MS2_VIEWER_TIME_RESOLUTION,threeDViewerMs2DefaultTimeResolution_);
    if (MS2_VIEWER_TIME_RESOLUTION.length()>longestKey) longestKey = MS2_VIEWER_TIME_RESOLUTION.length();
    if (threeDViewerMs2DefaultTimeResolution_.length()>longestValue) longestValue = threeDViewerMs2DefaultTimeResolution_.length(); 
    rowCount = createPropertyRow(sheet,rowCount,MS2_VIEWER_MZ_RESOLUTION,threeDViewerMs2DefaultMZResolution_);
    if (MS2_VIEWER_MZ_RESOLUTION.length()>longestKey) longestKey = MS2_VIEWER_MZ_RESOLUTION.length();
    if (threeDViewerMs2DefaultMZResolution_.length()>longestValue) longestValue = threeDViewerMs2DefaultMZResolution_.length();
    rowCount = createPropertyRow(sheet,rowCount,MAX_CHROM_AT_ONCE_MB,String.valueOf(maxFileSizeForChromTranslationAtOnceInMB_));
    if (MAX_CHROM_AT_ONCE_MB.length()>longestKey) longestKey = MAX_CHROM_AT_ONCE_MB.length();
    if (String.valueOf(maxFileSizeForChromTranslationAtOnceInMB_).length()>longestValue) longestValue = String.valueOf(maxFileSizeForChromTranslationAtOnceInMB_).length();   
    rowCount = createPropertyRow(sheet,rowCount,CHROM_MULT_FOR_INT,String.valueOf(chromMultiplicationFactorForInt_));
    if (CHROM_MULT_FOR_INT.length()>longestKey) longestKey = CHROM_MULT_FOR_INT.length();
    if (String.valueOf(chromMultiplicationFactorForInt_).length()>longestValue) longestValue = String.valueOf(chromMultiplicationFactorForInt_).length();   
    rowCount = createPropertyRow(sheet,rowCount,CHROM_RESOLUTION_LOWEST,String.valueOf(chromLowestResolution_));
    if (CHROM_RESOLUTION_LOWEST.length()>longestKey) longestKey = CHROM_RESOLUTION_LOWEST.length();
    if (String.valueOf(chromLowestResolution_).length()>longestValue) longestValue = String.valueOf(chromLowestResolution_).length();

      
    if (shotgun_>SHOTGUN_FALSE){
      rowCount = createPropertyRow(sheet,rowCount,SHOTGUN,String.valueOf(shotgun_));
      if (SHOTGUN.length()>longestKey) longestKey = SHOTGUN.length();
      if (String.valueOf(shotgun_).length()>longestValue) longestValue = String.valueOf(shotgun_).length();
      rowCount = createPropertyRow(sheet,rowCount,MZUNIT,mzUnit_);
      if (MZUNIT.length()>longestKey) longestKey = MZUNIT.length();
      if (mzUnit_.length()>longestValue) longestValue = mzUnit_.length();
    }
    if (shotgun_==SHOTGUN_TRUE){
      String shotgunProcessingString = null;
      if (shotgunProcessing_ == LipidomicsAnalyzer.SHOTGUN_TYPE_MEAN)
        shotgunProcessingString = SHOTGUN_PROCESSING_MEAN;
      else if (shotgunProcessing_ == LipidomicsAnalyzer.SHOTGUN_TYPE_MEDIAN)
        shotgunProcessingString = SHOTGUN_PROCESSING_MEDIAN;
      else if (shotgunProcessing_ == LipidomicsAnalyzer.SHOTGUN_TYPE_SUM)
        shotgunProcessingString = SHOTGUN_PROCESSING_SUM;
      rowCount = createPropertyRow(sheet,rowCount,SHOTGUN_PROCESSING,shotgunProcessingString);
      if (SHOTGUN_PROCESSING.length()>longestKey) longestKey = SHOTGUN_PROCESSING.length();
      if (shotgunProcessingString.length()>longestValue) longestValue = shotgunProcessingString.length();
      
      String shotgunZeroHandlingString = "false";
      if (shotgunRelIntCutoff_>0)
        shotgunZeroHandlingString = String.valueOf(shotgunRelIntCutoff_);
      else if (shotgunRelIntCutoff_==0)
        shotgunZeroHandlingString = "true";
      rowCount = createPropertyRow(sheet,rowCount,SHOTGUN_ZERO_HANDLING,shotgunZeroHandlingString);
      if (SHOTGUN_ZERO_HANDLING.length()>longestKey) longestKey = SHOTGUN_ZERO_HANDLING.length();
      if (shotgunZeroHandlingString.length()>longestValue) longestValue = shotgunZeroHandlingString.length();

    }else{
      rowCount = createPropertyRow(sheet,rowCount,CHROM_SMOOTH_RANGE,String.valueOf(chromSmoothRange_));
      if (CHROM_SMOOTH_RANGE.length()>longestKey) longestKey = CHROM_SMOOTH_RANGE.length();
      if (String.valueOf(chromSmoothRange_).length()>longestValue) longestValue = String.valueOf(chromSmoothRange_).length();
      rowCount = createPropertyRow(sheet,rowCount,CHROM_SMOOTH_REPEATS,String.valueOf(chromSmoothRepeats_));
      if (CHROM_SMOOTH_REPEATS.length()>longestKey) longestKey = CHROM_SMOOTH_REPEATS.length();
      if (String.valueOf(chromSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(chromSmoothRepeats_).length();
    }
    if (shotgun_==SHOTGUN_FALSE){
      rowCount = createPropertyRow(sheet,rowCount,USE_3D,String.valueOf(use3D_));
      if (USE_3D.length()>longestKey) longestKey = USE_3D.length();
      if (String.valueOf(use3D_).length()>longestValue) longestValue = String.valueOf(use3D_).length();
      rowCount = createPropertyRow(sheet,rowCount,ISOTOPE_CORRECTION,String.valueOf(isotopeCorrection_));
      if (ISOTOPE_CORRECTION.length()>longestKey) longestKey = ISOTOPE_CORRECTION.length();    
      if (String.valueOf(isotopeCorrection_).length()>longestValue) longestValue = String.valueOf(isotopeCorrection_).length();
      rowCount = createPropertyRow(sheet,rowCount,REMOVE_FROM_OTHER_ISOTOPE,String.valueOf(removeFromOtherIsotopes_));
      if (REMOVE_FROM_OTHER_ISOTOPE.length()>longestKey) longestKey = REMOVE_FROM_OTHER_ISOTOPE.length();
      if (String.valueOf(removeFromOtherIsotopes_).length()>longestValue) longestValue = String.valueOf(removeFromOtherIsotopes_).length();
      rowCount = createPropertyRow(sheet,rowCount,RESPECT_ISO_DISTRI,String.valueOf(respectIsotopicDistribution_));
      if (RESPECT_ISO_DISTRI.length()>longestKey) longestKey = RESPECT_ISO_DISTRI.length();
      if (String.valueOf(respectIsotopicDistribution_).length()>longestValue) longestValue = String.valueOf(respectIsotopicDistribution_).length();
      rowCount = createPropertyRow(sheet,rowCount,CHECK_CHAIN_LABEL_COMBINATION,String.valueOf(checkChainLabelCombination_));
      if (CHECK_CHAIN_LABEL_COMBINATION.length()>longestKey) longestKey = CHECK_CHAIN_LABEL_COMBINATION.length();
      if (String.valueOf(checkChainLabelCombination_).length()>longestValue) longestValue = String.valueOf(checkChainLabelCombination_).length();
      rowCount = createPropertyRow(sheet,rowCount,NOISE_CUTOFF,String.valueOf(useNoiseCutoff_));
      if (NOISE_CUTOFF.length()>longestKey) longestKey = NOISE_CUTOFF.length();
      if (String.valueOf(useNoiseCutoff_).length()>longestValue) longestValue = String.valueOf(useNoiseCutoff_).length();
      rowCount = createPropertyRow(sheet,rowCount,NOISE_DEVIATION,String.valueOf(noiseCutoffDeviationValue_));
      if (NOISE_DEVIATION.length()>longestKey) longestKey = NOISE_DEVIATION.length();
      if (String.valueOf(noiseCutoffDeviationValue_).length()>longestValue) longestValue = String.valueOf(noiseCutoffDeviationValue_).length();
      if (minimumRelativeIntensity_!=null) rowCount = createPropertyRow(sheet,rowCount,NOISE_MIN_INTENSITY,minimumRelativeIntensity_.toString());
      rowCount = createPropertyRow(sheet,rowCount,SCAN_STEP,String.valueOf(scanStep_));
      if (SCAN_STEP.length()>longestKey) longestKey = SCAN_STEP.length();
      if (String.valueOf(scanStep_).length()>longestValue) longestValue = String.valueOf(scanStep_).length();  
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_MZ_RANGE,String.valueOf(profileMzRange_*2));
      if (PROFILE_MZ_RANGE.length()>longestKey) longestKey = PROFILE_MZ_RANGE.length();
      if (String.valueOf(profileMzRange_*2).length()>longestValue) longestValue = String.valueOf(profileMzRange_*2).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_TIME_TOL,String.valueOf(profileTimeTolerance_));
      if (PROFILE_TIME_TOL.length()>longestKey) longestKey = PROFILE_TIME_TOL.length();
      if (String.valueOf(profileTimeTolerance_).length()>longestValue) longestValue = String.valueOf(profileTimeTolerance_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_INT_THRESHOLD,String.valueOf(profileIntThreshold_));
      if (PROFILE_INT_THRESHOLD.length()>longestKey) longestKey = PROFILE_INT_THRESHOLD.length();
      if (String.valueOf(profileIntThreshold_).length()>longestValue) longestValue = String.valueOf(profileIntThreshold_).length();   
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_BROADER_TIME_TOL,String.valueOf(broaderProfileTimeTolerance_));
      if (PROFILE_BROADER_TIME_TOL.length()>longestKey) longestKey = PROFILE_BROADER_TIME_TOL.length();
      if (String.valueOf(broaderProfileTimeTolerance_).length()>longestValue) longestValue = String.valueOf(broaderProfileTimeTolerance_).length();   
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_SMOOTH_RANGE,String.valueOf(profileSmoothRange_));
      if (PROFILE_SMOOTH_RANGE.length()>longestKey) longestKey = PROFILE_SMOOTH_RANGE.length();
      if (String.valueOf(profileSmoothRange_).length()>longestValue) longestValue = String.valueOf(profileSmoothRange_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_SMOOTH_REPEATS,String.valueOf(profileSmoothRepeats_));
      if (PROFILE_SMOOTH_REPEATS.length()>longestKey) longestKey = PROFILE_SMOOTH_REPEATS.length();
      if (String.valueOf(profileSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(profileSmoothRepeats_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_MEAN_SMOOTH_REPEATS,String.valueOf(profileMeanSmoothRepeats_));
      if (PROFILE_MEAN_SMOOTH_REPEATS.length()>longestKey) longestKey = PROFILE_MEAN_SMOOTH_REPEATS.length();
      if (String.valueOf(profileMeanSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(profileMeanSmoothRepeats_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_MZ_MIN_RANGE,String.valueOf(profileMzMinRange_));
      if (PROFILE_MZ_MIN_RANGE.length()>longestKey) longestKey = PROFILE_MZ_MIN_RANGE.length();
      if (String.valueOf(profileMzMinRange_).length()>longestValue) longestValue = String.valueOf(profileMzMinRange_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_STEEPNESS_CHANGE_1,String.valueOf(profileSteepnessChange1_));
      if (PROFILE_STEEPNESS_CHANGE_1.length()>longestKey) longestKey = PROFILE_STEEPNESS_CHANGE_1.length();
      if (String.valueOf(profileSteepnessChange1_).length()>longestValue) longestValue = String.valueOf(profileSteepnessChange1_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_STEEPNESS_CHANGE_2,String.valueOf(profileSteepnessChange2_));
      if (PROFILE_STEEPNESS_CHANGE_2.length()>longestKey) longestKey = PROFILE_STEEPNESS_CHANGE_2.length();
      if (String.valueOf(profileSteepnessChange2_).length()>longestValue) longestValue = String.valueOf(profileSteepnessChange2_).length();   
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_INT_CUTOFF_1,String.valueOf(profileIntensityCutoff1_));
      if (PROFILE_INT_CUTOFF_1.length()>longestKey) longestKey = PROFILE_INT_CUTOFF_1.length();
      if (String.valueOf(profileIntensityCutoff1_).length()>longestValue) longestValue = String.valueOf(profileIntensityCutoff1_).length();   
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_INT_CUTOFF_2,String.valueOf(profileIntensityCutoff2_));
      if (PROFILE_INT_CUTOFF_2.length()>longestKey) longestKey = PROFILE_INT_CUTOFF_2.length();
      if (String.valueOf(profileIntensityCutoff2_).length()>longestValue) longestValue = String.valueOf(profileIntensityCutoff2_).length();      
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_GENERAL_INT_CUTOFF,String.valueOf(profileGeneralIntCutoff_));
      if (PROFILE_GENERAL_INT_CUTOFF.length()>longestKey) longestKey = PROFILE_GENERAL_INT_CUTOFF.length();
      if (String.valueOf(profileGeneralIntCutoff_).length()>longestValue) longestValue = String.valueOf(profileGeneralIntCutoff_).length();      
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_PEAK_ACCEPTANCE_RANGE,String.valueOf(profilePeakAcceptanceRange_));
      if (PROFILE_PEAK_ACCEPTANCE_RANGE.length()>longestKey) longestKey = PROFILE_PEAK_ACCEPTANCE_RANGE.length();
      if (String.valueOf(profilePeakAcceptanceRange_).length()>longestValue) longestValue = String.valueOf(profilePeakAcceptanceRange_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_SMOOTHING_CORRECTION,String.valueOf(profileSmoothingCorrection_));
      if (PROFILE_SMOOTHING_CORRECTION.length()>longestKey) longestKey = PROFILE_SMOOTHING_CORRECTION.length();
      if (String.valueOf(profileSmoothingCorrection_).length()>longestValue) longestValue = String.valueOf(profileSmoothingCorrection_).length();
      rowCount = createPropertyRow(sheet,rowCount,PROFILE_MZ_MAX_RANGE,String.valueOf(profileMaxRange_));
      if (PROFILE_MZ_MAX_RANGE.length()>longestKey) longestKey = PROFILE_MZ_MAX_RANGE.length();
      if (String.valueOf(profileMaxRange_).length()>longestValue) longestValue = String.valueOf(profileMaxRange_).length();  
      rowCount = createPropertyRow(sheet,rowCount,SMALL_CHROM_MZ_RANGE,String.valueOf(smallChromMzRange_));
      if (SMALL_CHROM_MZ_RANGE.length()>longestKey) longestKey = SMALL_CHROM_MZ_RANGE.length();
      if (String.valueOf(smallChromMzRange_).length()>longestValue) longestValue = String.valueOf(smallChromMzRange_).length();  
      rowCount = createPropertyRow(sheet,rowCount,SMALL_CHROM_SMOOTH_REPEATS,String.valueOf(smallChromSmoothRepeats_));
      if (SMALL_CHROM_SMOOTH_REPEATS.length()>longestKey) longestKey = SMALL_CHROM_SMOOTH_REPEATS.length();
      if (String.valueOf(smallChromSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(smallChromSmoothRepeats_).length();      
      rowCount = createPropertyRow(sheet,rowCount,SMALL_CHROM_MEAN_SMOOTH_REPEATS,String.valueOf(smallChromMeanSmoothRepeats_));
      if (SMALL_CHROM_MEAN_SMOOTH_REPEATS.length()>longestKey) longestKey = SMALL_CHROM_MEAN_SMOOTH_REPEATS.length();
      if (String.valueOf(smallChromMeanSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(smallChromMeanSmoothRepeats_).length();
      rowCount = createPropertyRow(sheet,rowCount,SMALL_CHROM_SMOOTH_RANGE,String.valueOf(smallChromSmoothRange_));
      if (SMALL_CHROM_SMOOTH_RANGE.length()>longestKey) longestKey = SMALL_CHROM_SMOOTH_RANGE.length();
      if (String.valueOf(smallChromSmoothRange_).length()>longestValue) longestValue = String.valueOf(smallChromSmoothRange_).length();
      rowCount = createPropertyRow(sheet,rowCount,SMALL_CHROM_INT_CUTOFF,String.valueOf(smallChromIntensityCutoff_));
      if (SMALL_CHROM_INT_CUTOFF.length()>longestKey) longestKey = SMALL_CHROM_INT_CUTOFF.length();
      if (String.valueOf(smallChromIntensityCutoff_).length()>longestValue) longestValue = String.valueOf(smallChromIntensityCutoff_).length();
      rowCount = createPropertyRow(sheet,rowCount,BROAD_CHROM_SMOOTH_REPEATS,String.valueOf(broadChromSmoothRepeats_));
      if (BROAD_CHROM_SMOOTH_REPEATS.length()>longestKey) longestKey = BROAD_CHROM_SMOOTH_REPEATS.length();
      if (String.valueOf(broadChromSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(broadChromSmoothRepeats_).length();
      rowCount = createPropertyRow(sheet,rowCount,BROAD_CHROM_MEAN_SMOOTH_REPEATS,String.valueOf(broadChromMeanSmoothRepeats_));
      if (BROAD_CHROM_MEAN_SMOOTH_REPEATS.length()>longestKey) longestKey = BROAD_CHROM_MEAN_SMOOTH_REPEATS.length();
      if (String.valueOf(broadChromMeanSmoothRepeats_).length()>longestValue) longestValue = String.valueOf(broadChromMeanSmoothRepeats_).length();
      rowCount = createPropertyRow(sheet,rowCount,BROAD_CHROM_SMOOTH_RANGE,String.valueOf(broadChromSmoothRange_));
      if (BROAD_CHROM_SMOOTH_RANGE.length()>longestKey) longestKey = BROAD_CHROM_SMOOTH_RANGE.length();
      if (String.valueOf(broadChromSmoothRange_).length()>longestValue) longestValue = String.valueOf(broadChromSmoothRange_).length();
      rowCount = createPropertyRow(sheet,rowCount,BROAD_CHROM_INT_CUTOFF,String.valueOf(broadChromIntensityCutoff_));
      if (BROAD_CHROM_INT_CUTOFF.length()>longestKey) longestKey = BROAD_CHROM_INT_CUTOFF.length();
      if (String.valueOf(broadChromIntensityCutoff_).length()>longestValue) longestValue = String.valueOf(broadChromIntensityCutoff_).length();
      rowCount = createPropertyRow(sheet,rowCount,BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL,String.valueOf(broadChromSteepnessChangeNoSmall_));
      if (BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL.length()>longestKey) longestKey = BROAD_CHROM_STEEPNESS_CHANGE_NO_SMALL.length();
      if (String.valueOf(broadChromSteepnessChangeNoSmall_).length()>longestValue) longestValue = String.valueOf(broadChromSteepnessChangeNoSmall_).length();  
      rowCount = createPropertyRow(sheet,rowCount,BROAD_CHROM_INT_CUTOFF_NO_SMALL,String.valueOf(broadChromIntensityCutoffNoSmall_));
      if (BROAD_CHROM_INT_CUTOFF_NO_SMALL.length()>longestKey) longestKey = BROAD_CHROM_INT_CUTOFF_NO_SMALL.length();
      if (String.valueOf(broadChromIntensityCutoffNoSmall_).length()>longestValue) longestValue = String.valueOf(broadChromIntensityCutoffNoSmall_).length();  
      rowCount = createPropertyRow(sheet,rowCount,FINAL_PROBE_TIME_COMP_TOL,String.valueOf(finalProbeTimeCompTolerance_));
      if (FINAL_PROBE_TIME_COMP_TOL.length()>longestKey) longestKey = FINAL_PROBE_TIME_COMP_TOL.length();
      if (String.valueOf(finalProbeTimeCompTolerance_).length()>longestValue) longestValue = String.valueOf(finalProbeTimeCompTolerance_).length();   
      rowCount = createPropertyRow(sheet,rowCount,FINAL_PROBE_MZ_COMP_TOL,String.valueOf(finalProbeMzCompTolerance_));
      if (FINAL_PROBE_MZ_COMP_TOL.length()>longestKey) longestKey = FINAL_PROBE_MZ_COMP_TOL.length();
      if (String.valueOf(finalProbeMzCompTolerance_).length()>longestValue) longestValue = String.valueOf(finalProbeMzCompTolerance_).length();    
      rowCount = createPropertyRow(sheet,rowCount,OVERLAP_DIST_DEV_FACTOR,String.valueOf(overlapDistanceDeviationFactor_));
      if (OVERLAP_DIST_DEV_FACTOR.length()>longestKey) longestKey = OVERLAP_DIST_DEV_FACTOR.length();
      if (String.valueOf(overlapDistanceDeviationFactor_).length()>longestValue) longestValue = String.valueOf(overlapDistanceDeviationFactor_).length();
      rowCount = createPropertyRow(sheet,rowCount,OVERLAP_INT_THRESHOLD,String.valueOf(overlapPossibleIntensityThreshold_));
      if (OVERLAP_INT_THRESHOLD.length()>longestKey) longestKey = OVERLAP_INT_THRESHOLD.length();
      if (String.valueOf(overlapPossibleIntensityThreshold_).length()>longestValue) longestValue = String.valueOf(overlapPossibleIntensityThreshold_).length();
      rowCount = createPropertyRow(sheet,rowCount,OVERLAP_INT_SURE_THRESHOLD,String.valueOf(overlapSureIntensityThreshold_));
      if (OVERLAP_INT_SURE_THRESHOLD.length()>longestKey) longestKey = OVERLAP_INT_SURE_THRESHOLD.length();
      if (String.valueOf(overlapSureIntensityThreshold_).length()>longestValue) longestValue = String.valueOf(overlapSureIntensityThreshold_).length();
      rowCount = createPropertyRow(sheet,rowCount,OVERLAP_PEAK_DIST_DIVISOR,String.valueOf(overlapPeakDistanceDivisor_));
      if (OVERLAP_PEAK_DIST_DIVISOR.length()>longestKey) longestKey = OVERLAP_PEAK_DIST_DIVISOR.length();
      if (String.valueOf(overlapPeakDistanceDivisor_).length()>longestValue) longestValue = String.valueOf(overlapPeakDistanceDivisor_).length();
      rowCount = createPropertyRow(sheet,rowCount,OVERLAP_FULL_DIST_DIVISOR,String.valueOf(overlapFullDistanceDivisor_));
      if (OVERLAP_FULL_DIST_DIVISOR.length()>longestKey) longestKey = OVERLAP_FULL_DIST_DIVISOR.length();
      if (String.valueOf(overlapFullDistanceDivisor_).length()>longestValue) longestValue = String.valueOf(overlapFullDistanceDivisor_).length();
    }
    if (shotgun_!=SHOTGUN_TRUE){
      rowCount = createPropertyRow(sheet,rowCount,PEAK_DISCARD_AREA_FACTOR,String.valueOf(peakDiscardingAreaFactor_));
      if (PEAK_DISCARD_AREA_FACTOR.length()>longestKey) longestKey = PEAK_DISCARD_AREA_FACTOR.length();
      if (String.valueOf(peakDiscardingAreaFactor_).length()>longestValue) longestValue = String.valueOf(peakDiscardingAreaFactor_).length();
    }
    if (shotgun_==SHOTGUN_FALSE){
      rowCount = createPropertyRow(sheet,rowCount,ISO_IN_BETWEEN_TIME,String.valueOf(isotopeInBetweenTime_));
      if (ISO_IN_BETWEEN_TIME.length()>longestKey) longestKey = ISO_IN_BETWEEN_TIME.length();
      if (String.valueOf(isotopeInBetweenTime_).length()>longestValue) longestValue = String.valueOf(isotopeInBetweenTime_).length();
      rowCount = createPropertyRow(sheet,rowCount,ISO_IN_BETWEEN_AREA_FACTOR,String.valueOf(isoInBetweenAreaFactor_));
      if (ISO_IN_BETWEEN_AREA_FACTOR.length()>longestKey) longestKey = ISO_IN_BETWEEN_AREA_FACTOR.length();
      if (String.valueOf(isoInBetweenAreaFactor_).length()>longestValue) longestValue = String.valueOf(isoInBetweenAreaFactor_).length();
      rowCount = createPropertyRow(sheet,rowCount,ISO_NEAR_NORMAL_PROBE_TIME,String.valueOf(isoNearNormalProbeTime_));
      if (ISO_NEAR_NORMAL_PROBE_TIME.length()>longestKey) longestKey = ISO_NEAR_NORMAL_PROBE_TIME.length();
      if (String.valueOf(isoNearNormalProbeTime_).length()>longestValue) longestValue = String.valueOf(isoNearNormalProbeTime_).length();
    }
    if (shotgun_!=SHOTGUN_TRUE){
      rowCount = createPropertyRow(sheet,rowCount,RELATIVE_AREA_CUTOFF,String.valueOf(relativeAreaCutoff_));
      if (RELATIVE_AREA_CUTOFF.length()>longestKey) longestKey = RELATIVE_AREA_CUTOFF.length();
      if (String.valueOf(relativeAreaCutoff_).length()>longestValue) longestValue = String.valueOf(relativeAreaCutoff_).length();   
      rowCount = createPropertyRow(sheet,rowCount,RELATIVE_AREA_FAR_CUTOFF,String.valueOf(relativeFarAreaCutoff_));
      if (RELATIVE_AREA_FAR_CUTOFF.length()>longestKey) longestKey = RELATIVE_AREA_FAR_CUTOFF.length();
      if (String.valueOf(relativeFarAreaCutoff_).length()>longestValue) longestValue = String.valueOf(relativeFarAreaCutoff_).length();  
      rowCount = createPropertyRow(sheet,rowCount,RELATIVE_AREA_FAR_TIME_SPACE,String.valueOf(relativeFarAreaTimeSpace_));
      if (RELATIVE_AREA_FAR_TIME_SPACE.length()>longestKey) longestKey = RELATIVE_AREA_FAR_TIME_SPACE.length();
      if (String.valueOf(relativeFarAreaTimeSpace_).length()>longestValue) longestValue = String.valueOf(relativeFarAreaTimeSpace_).length();
    }
    if (shotgun_==SHOTGUN_FALSE){
      rowCount = createPropertyRow(sheet,rowCount,RELATIVE_ISO_INBETWEEN_CUTOFF,String.valueOf(relativeIsoInBetweenCutoff_));
      if (RELATIVE_ISO_INBETWEEN_CUTOFF.length()>longestKey) longestKey = RELATIVE_ISO_INBETWEEN_CUTOFF.length();
      if (String.valueOf(relativeIsoInBetweenCutoff_).length()>longestValue) longestValue = String.valueOf(relativeIsoInBetweenCutoff_).length();
      rowCount = createPropertyRow(sheet,rowCount,ISO_IN_BETWEEN_TIME_MAX,String.valueOf(isoInBetweenMaxTimeDistance_));
      if (ISO_IN_BETWEEN_TIME_MAX.length()>longestKey) longestKey = ISO_IN_BETWEEN_TIME_MAX.length();
      if (String.valueOf(isoInBetweenMaxTimeDistance_).length()>longestValue) longestValue = String.valueOf(isoInBetweenMaxTimeDistance_).length();
      rowCount = createPropertyRow(sheet,rowCount,TWIN_PEAK_MZ_TOL,String.valueOf(twinPeakMzTolerance_));
      if (TWIN_PEAK_MZ_TOL.length()>longestKey) longestKey = TWIN_PEAK_MZ_TOL.length();
      if (String.valueOf(twinPeakMzTolerance_).length()>longestValue) longestValue = String.valueOf(twinPeakMzTolerance_).length(); 
      rowCount = createPropertyRow(sheet,rowCount,PEAK_CLOSE_TIME_TOL,String.valueOf(closePeakTimeTolerance_));
      if (PEAK_CLOSE_TIME_TOL.length()>longestKey) longestKey = PEAK_CLOSE_TIME_TOL.length();
      if (String.valueOf(closePeakTimeTolerance_).length()>longestValue) longestValue = String.valueOf(closePeakTimeTolerance_).length(); 
      rowCount = createPropertyRow(sheet,rowCount,TWIN_INBETWEEN_CUTOFF,String.valueOf(twinInBetweenCutoff_));
      if (TWIN_INBETWEEN_CUTOFF.length()>longestKey) longestKey = TWIN_INBETWEEN_CUTOFF.length();
      if (String.valueOf(twinInBetweenCutoff_).length()>longestValue) longestValue = String.valueOf(twinInBetweenCutoff_).length();
      rowCount = createPropertyRow(sheet,rowCount,UNION_INBETWEEN_CUTOFF,String.valueOf(unionInBetweenCutoff_));
      if (UNION_INBETWEEN_CUTOFF.length()>longestKey) longestKey = UNION_INBETWEEN_CUTOFF.length();
      if (String.valueOf(unionInBetweenCutoff_).length()>longestValue) longestValue = String.valueOf(unionInBetweenCutoff_).length();
      rowCount = createPropertyRow(sheet,rowCount,SPARSE_DATA,String.valueOf(sparseData_));
      if (SPARSE_DATA.length()>longestKey) longestKey = SPARSE_DATA.length();
      if (String.valueOf(sparseData_).length()>longestValue) longestValue = String.valueOf(sparseData_).length();
    }
    if (mzTabInstrumentName_!=null){
      String value = getStringFromMzTabParameter(mzTabInstrumentName_);
      if (MZTAB_INSTRUMENT.length()>longestKey) longestKey = MZTAB_INSTRUMENT.length();
      if (value.length()>longestValue) longestValue = value.length();
      rowCount = createPropertyRow(sheet,rowCount,MZTAB_INSTRUMENT,value);
    }
    if (mzTabInstrumentSource_!=null){
      String value = getStringFromMzTabParameter(mzTabInstrumentSource_);
      if (MZTAB_IONSOURCE.length()>longestKey) longestKey = MZTAB_IONSOURCE.length();
      if (value.length()>longestValue) longestValue = value.length();
      rowCount = createPropertyRow(sheet,rowCount,MZTAB_IONSOURCE,value);
    }
    if (mzTabInstrumentAnalyzer_!=null){
      String value = getStringFromMzTabParameter(mzTabInstrumentAnalyzer_);
      if (MZTAB_MSANALYZER.length()>longestKey) longestKey = MZTAB_MSANALYZER.length();
      if (value.length()>longestValue) longestValue = value.length();
      rowCount = createPropertyRow(sheet,rowCount,MZTAB_MSANALYZER,value);
    }
    if (mzTabInstrumentDetector_!=null){
      String value = getStringFromMzTabParameter(mzTabInstrumentDetector_);
      if (MZTAB_DETECTOR.length()>longestKey) longestKey = MZTAB_DETECTOR.length();
      if (value.length()>longestValue) longestValue = value.length();
      rowCount = createPropertyRow(sheet,rowCount,MZTAB_DETECTOR,value);
    }
    if (alexTargetlist_==true){
      rowCount = createPropertyRow(sheet,rowCount,ALEX_TARGETLIST,String.valueOf(alexTargetlist_));
      if (ALEX_TARGETLIST.length()>longestKey) longestKey = ALEX_TARGETLIST.length();
      if (String.valueOf(alexTargetlist_).length()>longestValue) longestValue = String.valueOf(alexTargetlist_).length();
    }
    if (useMsconvertForWaters_==true) {
      rowCount = createPropertyRow(sheet,rowCount,USE_MSCONVERT_FOR_WATERS,String.valueOf(useMsconvertForWaters_));
      if (USE_MSCONVERT_FOR_WATERS.length()>longestKey) longestKey = USE_MSCONVERT_FOR_WATERS.length();
      if (String.valueOf(useMsconvertForWaters_).length()>longestValue) longestValue = String.valueOf(useMsconvertForWaters_).length();
    }
    
    
    String key;
    String value;
    if (faHydroxyEncoding!=null) {
      for (Short oh : faHydroxyEncoding.getHydroxyNumbersInAscendingOrder()) {
        try {
          key = EXCEL_HYDROXY_FA_PREFIX+String.valueOf(oh);
          value = faHydroxyEncoding.getEncodedPrefix(oh);
          rowCount = createPropertyRow(sheet,rowCount,key,value);
        //this catch can never happen
        }catch (HydroxylationEncodingException e) {}
      }
    }
    if (lcbHydroxyEncoding!=null) {
      for (Short oh : lcbHydroxyEncoding.getHydroxyNumbersInAscendingOrder()) {
        try {
          key = EXCEL_HYDROXY_LCB_PREFIX+String.valueOf(oh);
          value = lcbHydroxyEncoding.getEncodedPrefix(oh);
          rowCount = createPropertyRow(sheet,rowCount,key,value);
        //this catch can never happen
        }catch (HydroxylationEncodingException e) {}
      }
    }

    
    int keyColumnWidth = (int)((LipidomicsConstants.EXCEL_KEY.length()*ExcelUtils.CHAR_MULT)*ExcelUtils.BOLD_MULT);
    if ((longestKey+1)*ExcelUtils.CHAR_MULT>keyColumnWidth) keyColumnWidth =  (longestKey+1)*ExcelUtils.CHAR_MULT;
    sheet.setColumnWidth(EXCEL_KEY_COLUMN,keyColumnWidth); 
    int valueColumnWidth = (int)((LipidomicsConstants.EXCEL_VALUE.length()*ExcelUtils.CHAR_MULT)*ExcelUtils.BOLD_MULT);
    if ((longestValue+1)*ExcelUtils.CHAR_MULT>valueColumnWidth) valueColumnWidth =  (longestValue+1)*ExcelUtils.CHAR_MULT;
    sheet.setColumnWidth(EXCEL_VALUE_COLUMN,valueColumnWidth);

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

  /**
   * creates one Excel Row containing key and value pair
   * @param sheet the Excel sheet to write the row
   * @param rowCount index for the next available row
   * @param key key to be written (of the key/value pair)
   * @param value value to be written (of the key/value pair)
   * @return next available Excel row index
   */
  private int createPropertyRow(Sheet sheet,int rowCount,String key,String value){
    Row row = sheet.createRow(rowCount);
    rowCount++;
    setKeyCellValue(row,key);
    setValueCellValue(row,value);
    return rowCount;
  }
  
  /**
   * writes the key cell of the key/value pair
   * @param row row where the key has to be written
   * @param value the key to be written
   */
  private void setKeyCellValue(Row row, String value){
    setCellValue(row, EXCEL_KEY_COLUMN, value);
  }

  /**
   * writes the value cell of the key/value pair
   * @param row row where the key has to be written
   * @param value the value to be written
   */
  private void setValueCellValue(Row row, String value){
    setCellValue(row, EXCEL_VALUE_COLUMN, value);
  }
  
  /**
   * writes the value in a certain column of the defined row
   * @param row row where the cell has to be written
   * @param column column where the cell has to be written
   * @param value value to be written in the cell
   */
  private void setCellValue(Row row, int column, String value){
    Cell cell = row.createCell(column);
    cell.setCellValue(value);
  }
  
  /**
   * creates a LipidomicsConstants object from an Excel sheet - this cannot be done in a separate
   * class, since the parameters to be set are private
   * @param sheet Excel sheet to read the parameters from
   * @throws SettingsException thrown when a settings combination is not possible
   * @return settings: [0] LipidomicsConstants object containing the parameters that were read; [1] FA hydroxylation encoding; [2] LCB hydroxylation encoding
   */
  public static Object[] readSettingsFromExcel(Sheet sheet) throws SettingsException{
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
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_KEY)) keyColumn = columnId;
          if (entry.equalsIgnoreCase(LipidomicsConstants.EXCEL_VALUE)) valueColumn = columnId;
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

  public String getRelativeMS1BasePeakCutoff()
  {
    return relativeMS1BasePeakCutoff_;
  }

  public void setRelativeMS1BasePeakCutoff(String relativeMS1BasePeakCutoff)
  {
    this.relativeMS1BasePeakCutoff_ = relativeMS1BasePeakCutoff;
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
    if (instance_ == null) LipidomicsConstants.getInstance();
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
