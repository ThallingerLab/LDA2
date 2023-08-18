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

package at.tugraz.genome.lda.quantification;

//import java.io.BufferedOutputStream;
//import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

//import javax.imageio.ImageIO;
//
//import org.jfree.chart.ChartFactory;
//import org.jfree.chart.JFreeChart;
//import org.jfree.chart.plot.PlotOrientation;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;








import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.Analyzer;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.QuantificationException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.FragmentCalculator;
import at.tugraz.genome.lda.msn.vos.FragmentRuleVO;
import at.tugraz.genome.lda.msn.vos.FragmentVO;
import at.tugraz.genome.lda.quantification.LipidomicsDefines;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.ChromaAnalyzer;
import at.tugraz.genome.maspectras.quantification.Probe3D;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsAnalyzer extends ChromaAnalyzer
{
  
  private boolean useCuda_;
  protected SavGolJNI sav_gol_jni_;
  
  public LipidomicsAnalyzer(String headerFilePath, String indexFilePath, String retentionTimeFilePath, String chromatogramFilePath, boolean useCuda)throws CgException{
    super();
    if (useCuda){
      sav_gol_jni_ = new SavGolJNI();
    }
    reader_ = new LipidomicsChromReader(headerFilePath,indexFilePath,retentionTimeFilePath,chromatogramFilePath,LipidomicsConstants.isSparseData(),
        LipidomicsConstants.getChromSmoothRange(), useCuda);
    m_chroma = new LipidomicsChromatogram[CgDefines.MaxCharge];
    LipidomicsChromReader lReader = (LipidomicsChromReader) reader_;
    if (useCuda){
      int max_length = lReader.getNumberOfScans();
      for(Integer key : lReader.getMsmsNrOfScans().keySet()){
        if (max_length < lReader.getMsmsNrOfScans().get(key)){
          max_length = lReader.getMsmsNrOfScans().get(key);
        }
      }
      sav_gol_jni_.initMalloc( max_length );
    }
    this.useCuda_ = useCuda;
    this.init();
  }
  
  private boolean useSameCgHashFor3D_;
//  public static float NEUTRON_MASS = 1.00866491597f;
  // here just a mean value for the mass difference of one isotope is
  // taken, because an isotope has never the real mass difference of one
  // neutron; the reason is that neutrons in bound form are not so heavy
  // the mass difference in biological elements range normally from 
  // 0.999Da-1.007Da, however, the ones with lower masses than 1Da do
  // not occur so often!
  // This is no settable in the properties file
  //public static float NEUTRON_MASS = 1.0025f;
  
  private static int NO_ISOTOPE_PRESENT = 0;
  private static int POSSIBLE_ISOTOPE_OVERLAP = 1;
  private static int ISOTOPE_PRESENT = 2;
  
//  public boolean verbose_;
//  private Vector<CgProbe> sameCgHash_;
//  private Vector<Probe3D> same3DProbes_;


  /** the following parameters are necessary for lipidomics;
   *  since MASPECTRAS uses this class as well, I cannot simply 
   *  take them from the LipidomicsConstants, but I have to initialize
   *  them
   */
  private float generalBasePeakCutoff_;
  private float coarseChromMzTolerance_;
  private float chromSmoothRange_= 0.5f;
  private int chromSmoothRepeats_ = 10;
  private boolean removeIfOtherIsotopePresent_ = true; 
  private boolean useNoiseCutoff_;
  private float noiseCutoffDeviationValue_;
  private Float minimumRelativeIntensity_;
  private int scanStep_;
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
  private float broadChromSmoothRange_;
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
  @SuppressWarnings("unused")
  private float twinPeakMzTolerance_;
  private int closePeakTimeTolerance_;
  private float twinInBetweenCutoff_;
  private float unionInBetweenCutoff_;

  // this is for caching the chromatograms to speed up the process
  private Vector<CgProbe> chromProbes_;
  private Hashtable<Integer,CgProbe> fromProbeToProfile_;
  private Vector<CgProbe> profileProbes_;
  private Hashtable<Integer,LipidomicsChromatogram> fromProfileToSmallChrom_;
  private Hashtable<Integer,LipidomicsChromatogram> fromProfileToBroadChrom_;
  private Hashtable<Float,CgProbe> prevIsoProbes_;
  
  private float highestArea_ = 0f;
  private float highestIntensity_ = 0f;
  
  /** caching the chrom files - MS1 chrom files are used by isobar separation algorithm*/
  private  Hashtable<Integer,LipidomicsChromatogram> chromHash_;
  
//  private int calculateProfile_ = 0;
  
  private float lowerIntensityThreshold_ = 0f;
//  private Vector<CgProbe> smallProbes_;
//  private Vector<CgProbe> broadProbes_;
  
  private float msnMzTolerance_ = 0.2f;

  /** the shotgun intensity value is a mean of the scans*/
  public final static int SHOTGUN_TYPE_MEAN = 0;
  /** the shotgun intensity value is a median of the scans*/
  public final static int SHOTGUN_TYPE_MEDIAN = 1;
  /** the shotgun intensity value is the sum of the scan intensities*/
  public final static int SHOTGUN_TYPE_SUM = 2;
  /** which intensity type shall be used for shotgun data*/;
  private int shotgunType_;

  private void init(){
    useSameCgHashFor3D_ = false;
//    sameCgHash_ = new Vector<CgProbe>();
//    same3DProbes_ = new Vector<Probe3D>();
    this.initCacheHashes();
    this.generalBasePeakCutoff_ = 0.0f;
    this.coarseChromMzTolerance_ = 0.02f;
    this.removeIfOtherIsotopePresent_ = true;
    this.useNoiseCutoff_ = false;
    this.noiseCutoffDeviationValue_ = 2f;
    this.scanStep_ = 1;
    this.chromSmoothRange_= 0.5f;
    this.chromSmoothRepeats_ = 10;
    this.profileMzRange_ = 0.1f;
    this.profileTimeTolerance_ =  2f;
    this.profileIntThreshold_ = 5f;
    this.broaderProfileTimeTolerance_ = 3f;
    this.profileSmoothRange_ = 0.02f;
    this.profileSmoothRepeats_ = 5;
    this.profileMeanSmoothRepeats_ = 0;
    this.profileMzMinRange_ = 0.0f;
    this.profileSteepnessChange1_ = 1.5f;
    this.profileSteepnessChange2_ = 1.8f;
    this.profileIntensityCutoff1_ = 0.15f;
    this.profileIntensityCutoff2_ = 0.2f;
    this.profileGeneralIntCutoff_ = 0.03f;
    this.profilePeakAcceptanceRange_ = 0.034f;
    this.profileSmoothingCorrection_ = 0.0f;
    this.profileMaxRange_ = 0.02f;
    this.smallChromMzRange_ = 0.002f;
    this.smallChromSmoothRepeats_ = 15;
    this.smallChromMeanSmoothRepeats_ = 0;
    this.smallChromSmoothRange_ = 5f;
    this.smallChromIntensityCutoff_ = 0.03f;
    this.broadChromSmoothRepeats_ = 10;
    this.broadChromMeanSmoothRepeats_ = 0;
    this.broadChromSmoothRange_ = 5f;
    this.broadChromIntensityCutoff_= 0f;
    this.broadChromSteepnessChangeNoSmall_ = 1.33f;
    this.broadChromIntensityCutoffNoSmall_ = 0.05f;
    this.finalProbeTimeCompTolerance_ = 0.1f;
    this.finalProbeMzCompTolerance_ = 0.0005f;
    this.overlapDistanceDeviationFactor_ = 1.5f;
    this.overlapPossibleIntensityThreshold_ = 0.15f;
    this.overlapSureIntensityThreshold_ = 0.7f;
    this.overlapPeakDistanceDivisor_ = 3f;
    this.overlapFullDistanceDivisor_ = 6f;
    this.peakDiscardingAreaFactor_= 1000;
    this.isotopeInBetweenTime_ = 30;
    this.isoInBetweenAreaFactor_ = 3f;
    this.isoNearNormalProbeTime_ = 30;
    this.relativeAreaCutoff_ = 0.01f;
    this.relativeFarAreaCutoff_ = 0.1f;
    this.relativeFarAreaTimeSpace_ = 30;
    this.relativeIsoInBetweenCutoff_ = 0.5f;
    this.twinPeakMzTolerance_ = 0.002f;
    this.closePeakTimeTolerance_ = 10;
    this.isoInBetweenMaxTimeDistance_ = 300;
    this.twinInBetweenCutoff_ = 0.95f;
    this.unionInBetweenCutoff_ = 0.8f;
    
    this.msnMzTolerance_ = 0.2f;
  }
  
  public void set3DParameters(float chromSmoothRange, int chromSmoothRepeats, boolean removeIfOtherIsotopePresent,
      boolean useNoiseCutoff, float noiseCutoffDeviationValue, Float minimumRelativeIntensity, int scanStep, float profileMzRange, float profileTimeTolerance,
      float profileIntThreshold, float broaderProfileTimeTolerance,
      float profileSmoothRange, int profileSmoothRepeats, int profileMeanSmoothRepeats, float profileMzMinRange, float profileSteepnessChange1, float profileSteepnessChange2,
      float profileIntensityCutoff1, float profileIntensityCutoff2, float profileGeneralIntCutoff, float profilePeakAcceptanceRange, float profileSmoothingCorrection,
      float profileMaxRange, float smallChromMzRange,int smallChromSmoothRepeats, int smallChromMeanSmoothRepeats, float smallChromSmoothRange, float smallChromIntensityCutoff,
      int broadChromSmoothRepeats, int broadChromMeanSmoothRepeats, float broadChromSmoothRange, float broadChromIntensityCutoff, float broadChromSteepnessChangeNoSmall, 
      float broadChromIntensityCutoffNoSmall, float finalProbeTimeCompTolerance, float finalProbeMzCompTolerance,
      float overlapDistanceDeviationFactor, float overlapPossibleIntensityThreshold, float overlapSureIntensityThreshold,
      float overlapPeakDistanceDivisor, float overlapFullDistanceDivisor, int peakDiscardingAreaFactor,
      int isotopeInBetweenTime,  float isoInBetweenAreaFactor, int isoInBetweenMaxTimeDistance, int isoNearNormalProbeTime,
      float relativeAreaCutoff, float relativeFarAreaCutoff, int relativeFarAreaTimeSpace,
      float relativeIsoInBetweenCutoff, float twinPeakMzTolerance, int closePeakTimeTolerance, float twinInBetweenCutoff, float unionInBetweenCutoff,
      float msnMzTolerance){
    this.chromSmoothRange_= chromSmoothRange;
    this.removeIfOtherIsotopePresent_ = removeIfOtherIsotopePresent;
    this.chromSmoothRepeats_ = chromSmoothRepeats;
    this.useNoiseCutoff_ = useNoiseCutoff;
    this.noiseCutoffDeviationValue_ = noiseCutoffDeviationValue;
    this.minimumRelativeIntensity_ = minimumRelativeIntensity;
    this.scanStep_ = scanStep;
    this.profileMzRange_ = profileMzRange;
    this.profileTimeTolerance_ = profileTimeTolerance;
    this.profileIntThreshold_ = profileIntThreshold;
    this.broaderProfileTimeTolerance_ = broaderProfileTimeTolerance;
    this.profileSmoothRange_ = profileSmoothRange;
    this.profileSmoothRepeats_ = profileSmoothRepeats;
    this.profileMeanSmoothRepeats_ = profileMeanSmoothRepeats;
    this.profileMzMinRange_ = profileMzMinRange;
    this.profileSteepnessChange1_ = profileSteepnessChange1;
    this.profileSteepnessChange2_ = profileSteepnessChange2;
    this.profileIntensityCutoff1_ = profileIntensityCutoff1;
    this.profileIntensityCutoff2_ = profileIntensityCutoff2;
    this.profileGeneralIntCutoff_ = profileGeneralIntCutoff;
    this.profilePeakAcceptanceRange_ = profilePeakAcceptanceRange;
    this.profileSmoothingCorrection_ = profileSmoothingCorrection;
    this.profileMaxRange_ = profileMaxRange;
    this.smallChromMzRange_ = smallChromMzRange;
    this.smallChromSmoothRepeats_ = smallChromSmoothRepeats;
    this.smallChromMeanSmoothRepeats_ = smallChromMeanSmoothRepeats;
    this.smallChromSmoothRange_ = smallChromSmoothRange;
    this.smallChromIntensityCutoff_ = smallChromIntensityCutoff;
    this.broadChromSmoothRepeats_ = broadChromSmoothRepeats;
    this.broadChromMeanSmoothRepeats_ = broadChromMeanSmoothRepeats;
    this.broadChromSmoothRange_ = broadChromSmoothRange;
    this.broadChromIntensityCutoff_ = broadChromIntensityCutoff;
    this.broadChromSteepnessChangeNoSmall_ = broadChromSteepnessChangeNoSmall;
    this.broadChromIntensityCutoffNoSmall_ = broadChromIntensityCutoffNoSmall;
    this.finalProbeTimeCompTolerance_ = finalProbeTimeCompTolerance;
    this.finalProbeMzCompTolerance_ = finalProbeMzCompTolerance;
    this.overlapDistanceDeviationFactor_ = overlapDistanceDeviationFactor;
    this.overlapPossibleIntensityThreshold_ = overlapPossibleIntensityThreshold;
    this.overlapSureIntensityThreshold_ = overlapSureIntensityThreshold;
    this.overlapPeakDistanceDivisor_ = overlapPeakDistanceDivisor;
    this.overlapFullDistanceDivisor_ = overlapFullDistanceDivisor;
    this.peakDiscardingAreaFactor_ = peakDiscardingAreaFactor;
    this.isotopeInBetweenTime_ = isotopeInBetweenTime;
    this.isoInBetweenAreaFactor_ = isoInBetweenAreaFactor;
    this.isoInBetweenMaxTimeDistance_ = isoInBetweenMaxTimeDistance;
    this.isoNearNormalProbeTime_ = isoNearNormalProbeTime;
    this.relativeAreaCutoff_ = relativeAreaCutoff;
    this.relativeFarAreaCutoff_ = relativeFarAreaCutoff;
    this.relativeFarAreaTimeSpace_ = relativeFarAreaTimeSpace;
    this.relativeIsoInBetweenCutoff_ = relativeIsoInBetweenCutoff;
    this.twinPeakMzTolerance_ = twinPeakMzTolerance;
    this.closePeakTimeTolerance_ = closePeakTimeTolerance;
    this.twinInBetweenCutoff_ = twinInBetweenCutoff;
    this.unionInBetweenCutoff_ = unionInBetweenCutoff;
    this.msnMzTolerance_ = msnMzTolerance;
  }
  
  /**
   * sets parameters necessary for PRM procession
   * @param chromSmoothRange the range in seconds for smoothing the initial chromatogram
   * @param chromSmoothRepeats the smooth repeats for smoothing the initial chromatogram
   * @param peakDiscardingAreaFactor a factor for discarding small peaks immediately after detection, e.g. 100 corresponds to removal of peaks that are smaller than 1% of the strongest peak in the chromatogram
   * @param relativeAreaCutoff a factor for discarding peaks after all checks were performed, e.g. 0.01 corresponds to removal of peaks that are smaller than 1% of the strongest peak in the chromatogram
   * @param relativeFarAreaCutoff a factor for discarding peaks that are farer away from the strongest peak after all checks were performed, e.g. 0.1 corresponds to removal of peaks that are smaller than 10% of the strongest peak in the chromatogram
   * @param relativeFarAreaTimeSpace the time space in seconds to decide whether a peak is farer away from the strongest peak or not
   * @param msnMzTolerance the m/z tolerance in MSn spectra
   */
  public void setPrmParameters(float chromSmoothRange, int chromSmoothRepeats, int peakDiscardingAreaFactor, float relativeAreaCutoff,
      float relativeFarAreaCutoff, int relativeFarAreaTimeSpace, float msnMzTolerance){
    this.chromSmoothRange_= chromSmoothRange;
    this.chromSmoothRepeats_ = chromSmoothRepeats;
    this.peakDiscardingAreaFactor_ = peakDiscardingAreaFactor;
    this.relativeAreaCutoff_ = relativeAreaCutoff;
    this.relativeFarAreaCutoff_ = relativeFarAreaCutoff;
    this.relativeFarAreaTimeSpace_ = relativeFarAreaTimeSpace;
    this.msnMzTolerance_ = msnMzTolerance;
  }
  
  
  
  /**
   * sets the parameters necessary for processing shotgun data
   * @param mzTolerance the mzTolerance
   */
  public void setShotgunParameters(int shotgunType){
    this.shotgunType_ = shotgunType;
  }
    

  /** The hashes are used to reduce the amount of chromatogram reads*/
  private void initCacheHashes(){
    this.chromProbes_ = new Vector<CgProbe>();
    this.fromProbeToProfile_ = new Hashtable<Integer,CgProbe>();
    this.profileProbes_ = new Vector<CgProbe>();
    this.fromProfileToSmallChrom_ = new Hashtable<Integer,LipidomicsChromatogram>();
    this.fromProfileToBroadChrom_ = new Hashtable<Integer,LipidomicsChromatogram>();
    lowerIntensityThreshold_ = 0f;
//    this.smallProbes_ = new Vector<CgProbe>();
//    this.broadProbes_ = new Vector<CgProbe>();    
  }
  
  private void releaseCacheHashes(){
    this.chromProbes_ = null;
    this.fromProbeToProfile_ = null;
    this.profileProbes_ = null;
    this.fromProfileToSmallChrom_ = null;
    this.fromProfileToBroadChrom_ = null;
    lowerIntensityThreshold_ = 0f;
  }

  protected CgChromatogram readAChromatogram(float mz, float lowerMzBand, float upperMzBand,int msLevel,float smoothRange, int smoothRepeats, float meanSmoothRange, int meanSmoothRepeats, float startTime, float stopTime) throws CgException{
//  long time = System.currentTimeMillis();
	CgChromatogram cx = reader_.readChromatogram(mz - lowerMzBand,mz + upperMzBand,msLevel);
    cx.Mz = mz;
    cx.LowerMzBand = lowerMzBand;
    cx.UpperMzBand = upperMzBand;
    if (startTime!=-1f || stopTime!=Float.MAX_VALUE){
      boolean foundStart = false;
      boolean foundStop = false;
      for (int i=0;i!=cx.ScanCount;i++){
        if (!foundStart && startTime<cx.Value[i][0]){
          cx.startSmoothScan_ = i-1;
          foundStart = true;
        }
        if (!foundStop && stopTime<cx.Value[i][0]){
          cx.stopSmoothScan_ = i+1;
          foundStop = true;
        }
      }
    }
    boolean copyRawData = true;
    if (meanSmoothRange>0){
      cx.smoothMean(meanSmoothRange, meanSmoothRepeats,copyRawData);
      copyRawData = false;
    }
//    System.out.println("Reading Time: "+(System.currentTimeMillis()-time));
//    time = System.currentTimeMillis();
    if (useCuda_){
      LipidomicsChromatogram lCx = new LipidomicsChromatogram(cx);
      lCx.Smooth(smoothRange, smoothRepeats, copyRawData, sav_gol_jni_);
    } else {
      cx.Smooth(smoothRange, smoothRepeats, copyRawData);
    }
    
//    System.out.println("Smoothing Time: "+(System.currentTimeMillis()-time));
    return cx;
  }
  
  public Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processByMzProbabsAndPossibleRetentionTime(float mz, int charge, float retentionTime, float prevTimeTolerance, float afterTimeTolerance, int timeType, Vector<Double>probabs, Vector<Double>possibleProbabs,int msLevel, boolean negative) throws CgException{
    return this.processByMzProbabsAndPossibleRetentionTime(mz, charge, retentionTime, null, prevTimeTolerance, afterTimeTolerance, timeType, probabs, possibleProbabs, msLevel, negative);
  }

  public Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processByMzProbabsAndPossibleRetentionTime(float mz, int charge, Vector<Float> retentionTimes, float prevTimeTolerance, float afterTimeTolerance, int timeType, Vector<Double>probabs, Vector<Double>possibleProbabs,int msLevel, boolean negative) throws CgException{
    return this.processByMzProbabsAndPossibleRetentionTime(mz, charge, -1f, retentionTimes, prevTimeTolerance, afterTimeTolerance, timeType, probabs, possibleProbabs, msLevel, negative);
  }
  
  /**
   * 
   * @param mz the m/z value of the analyte
   * @param charge the charge of the analyte
   * @param msLevel the MS-level
   * @param nrOfIsotopes how many isotopes shall be quantified (starting with 0)
   * @return the detected results
   * @throws CgException
   */
  public Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processShotgunData(float mz, int charge, int msLevel, int nrOfIsotopes) throws CgException{
    float mzTolerance = LipidomicsConstants.getCoarseChromMzTolerance(mz);
    CgProbe probe = calculateAShotgunIntensity(mz,mzTolerance,charge,msLevel);
    if (!(probe.Area>0f))
      return null;
    Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> result = new Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>>();
    Hashtable<Integer,Vector<CgProbe>> oneResult = new  Hashtable<Integer,Vector<CgProbe>>();
    Vector<CgProbe> probes = new Vector<CgProbe>();
    probes.add(probe);
    oneResult.put(0, probes);
    if (nrOfIsotopes>1){
      for (int i=1; i!=nrOfIsotopes; i++){
        probe = calculateAShotgunIntensity(mz+i*LipidomicsConstants.getNeutronMass()/(float)charge,mzTolerance,charge,msLevel);
        if (probe.Area==0f)
          break;
        else{
          probes = new Vector<CgProbe>();
          probes.add(probe);
          oneResult.put(i, probes);
        }
      }
    }
    result.put(0, oneResult);
    return result;
  }
  
  /**
   * currently the method works only for head group fragments
   * @param mz the ideal m/z value
   * @param charge the charge
   * @param msLevel the MS-level
   * @param className name of the analyte class
   * @param modName name of modification
   * @param analyteName name of analyte
   * @param formula formula of analyte
   * @param ohNumber nuber of hydroxylation sites
   * @return quantification results
   */
  public Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processPrmData(float mz, int charge, int msLevel, String className, String modName,
      String analyteName, String formula, int ohNumber, String oxState) {
    Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> finalResults = null;
    try{
      float mzTolerance = LipidomicsConstants.getCoarseChromMzTolerance(mz);
      this.prepareMSnSpectraCache(mz-mzTolerance, mz+mzTolerance,LipidomicsConstants.getMs2MinIntsForNoiseRemoval());

      //TODO: here, the analyteFormula without deducts is not implemented - has to be changed in future!
      FragmentCalculator fragCalc = new FragmentCalculator(null,className,modName,analyteName,formula,formula,mz,charge,ohNumber,oxState);
      Vector<FragmentVO> mandHeadFragments = fragCalc.getHeadFragments(ohNumber).get(true);
      Hashtable<String,Vector<CgProbe>> headPeaks = new Hashtable<String,Vector<CgProbe>>();
      for (FragmentVO frag : mandHeadFragments){
        if (frag.isMandatory()!=FragmentRuleVO.MANDATORY_QUANT)
          continue;
        
        LipidomicsChromatogram chrom = new LipidomicsChromatogram(this.readJustIntensitiesOfInterest((float)frag.getMass()-msnMzTolerance_,(float)frag.getMass()+msnMzTolerance_,
            0f,Float.MAX_VALUE,msLevel));
        chrom.Smooth(chromSmoothRange_, this.chromSmoothRepeats_);
        chrom.GetMaximumAndAverage();

//        printChromaToFile(chrom,"D:\\Alex\\PRM\\test.png");
//        printChromaToFile(chrom,"D:\\Alex\\PRM\\test_raw.png",1);
        
        //getting the peaks of one chromatogram
        int currentScan =0;
        Vector<CgProbe> probes = new Vector<CgProbe>();
        float highestArea = 0f;
        while (currentScan<chrom.Value.length){
          CgProbe standardProbe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.StandardValleyMethod,charge);
          int highestScan = LipidomicsAnalyzer.findIndexByTime(standardProbe.UpperValley,chrom);
          standardProbe.isotopeNumber = 0;
          if (standardProbe.AreaStatus==CgAreaStatus.OK){
            this.checkDuplicate(probes,standardProbe);
            if (standardProbe.AreaStatus == CgAreaStatus.OK){
              probes.add(standardProbe);
              if (standardProbe.Area>highestArea)
                highestArea = standardProbe.Area;
              if (highestScan>currentScan)
              currentScan = highestScan;
            }
          }
          currentScan++;
        }
        //filter out small intensities
        Vector<CgProbe> strongProbes = new Vector<CgProbe>();
        for (CgProbe probe : probes){
          if (probe.Area>(highestArea/peakDiscardingAreaFactor_))
            strongProbes.add(probe);
        }
        headPeaks.put(frag.getName(), strongProbes);
      }
      finalResults = new Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>>();
    
      //TODO: here the chain Peaks have to be extracted - for the chain peaks, I have to check that all chains of one combination are present
      Hashtable<String,Vector<CgProbe>> chainPeaks = new Hashtable<String,Vector<CgProbe>>();
    
      if (headPeaks.size()==0 && chainPeaks.size()==0)
        return finalResults;
          
      //only hits which are found for all fragments are allowed
      //this method checks them and generates a final CgProbe including the areas
      String fragmentName;
      Vector<CgProbe> probesToCheck = new Vector<CgProbe>();
      if (headPeaks.size()>0){
        fragmentName = headPeaks.keySet().iterator().next();
        probesToCheck = headPeaks.get(fragmentName);
      }else{
        fragmentName = chainPeaks.keySet().iterator().next();
        probesToCheck = chainPeaks.get(fragmentName);
      }
      
      Vector<CgProbe> allProbes = new  Vector<CgProbe>();
      CgProbe strongestProbe = null;
      float strongestArea = 0f;
      for (CgProbe probe : probesToCheck){
        Vector<CgProbe> probesToCombine = new  Vector<CgProbe>();
        probesToCombine.add(probe);
        boolean allFound = true;
        for (String fragName : headPeaks.keySet()){
          if (fragName.equalsIgnoreCase(fragmentName))
            continue;
          Vector<CgProbe> otherProbes = headPeaks.get(fragName);
          Vector<CgProbe> overlappingProbes = new Vector<CgProbe>();
          CgProbe probeToAdd = null;
          for (CgProbe other: otherProbes){
            //here not inner third is used, but the peak of the one peak has to be within 75% of the other peak
            if (isOneProbeInOtherInnerThird(probe, other, 4/3f, 8/3f)){
              overlappingProbes.add(other);
            }
          }
          // if there is none found, discard the peak
          if (overlappingProbes.size()==0){
            allFound = false;
            break;
            // if  there is more than one found, take the one closest to the peak center
          } else {
            float lowestDifference = Float.MAX_VALUE;
            float summitDifference;
            for (CgProbe other : overlappingProbes){
              summitDifference = positiveDifference(probe.Peak, other.Peak);
              if (summitDifference<lowestDifference){
                probeToAdd = other;
                lowestDifference = summitDifference;
              }
            }
            probesToCombine.add(probeToAdd);
          }
        }
        if (!allFound)
          continue;
        //now make a combined area and add it to the final results
        float totalArea = 0f;
        float lowerValleyTimesArea = 0f;
        float peakTimesArea = 0f;
        float upperValleyTimesArea = 0;
        float totalBackground = 0f;
        float totalAreaError = 0f;
        float highestIntensity = 0f;
        float lowerValley10TimesArea = 0f;
        float lowerValley50TimesArea = 0f;
        float upperValley10TimesArea = 0f;
        float upperValley50TimesArea = 0f;
        for (CgProbe aProbe : probesToCombine){
          totalArea += aProbe.Area;
          lowerValleyTimesArea += aProbe.LowerValley*aProbe.Area;
          peakTimesArea += aProbe.Peak*aProbe.Area;
          upperValleyTimesArea += aProbe.UpperValley*aProbe.Area;
          totalBackground += aProbe.Background;
          totalAreaError += aProbe.AreaError;
          if (aProbe.getHighestIntensity()>highestIntensity)
            highestIntensity = aProbe.getHighestIntensity();
          lowerValley10TimesArea += aProbe.getLowerValley10()*aProbe.Area;
          lowerValley50TimesArea += aProbe.getLowerValley50()*aProbe.Area;
          upperValley10TimesArea += aProbe.getUpperValley10()*aProbe.Area;
          upperValley50TimesArea += aProbe.getUpperValley50()*aProbe.Area;
        }
        CgProbe finalProbe = new CgProbe(0,1,1,formula);
        finalProbe.Mz = mz;
        finalProbe.LowerMzBand = mzTolerance;
        finalProbe.UpperMzBand = mzTolerance;
        finalProbe.isotopeNumber = 0;
        finalProbe.Area = totalArea;
        finalProbe.LowerValley = lowerValleyTimesArea/totalArea;
        finalProbe.Peak = peakTimesArea/totalArea;
        finalProbe.UpperValley = upperValleyTimesArea/totalArea;
        finalProbe.Background = totalBackground;
        finalProbe.AreaError = totalAreaError;
        finalProbe.setHighestIntensity(highestIntensity);
        finalProbe.setLowerValley10(lowerValley10TimesArea/totalArea);
        finalProbe.setLowerValley50(lowerValley50TimesArea/totalArea);
        finalProbe.setUpperValley10(upperValley10TimesArea/totalArea);
        finalProbe.setUpperValley50(upperValley50TimesArea/totalArea);
        finalProbe.AreaStatus = CgAreaStatus.OK;
        finalProbe.setApexIntensity(0f);
        allProbes.add(finalProbe);
        if (finalProbe.Area>strongestArea){
          strongestProbe = finalProbe;
          strongestArea = finalProbe.Area;
        }
      }
      int hitNumber = 0;
      for (CgProbe aProbe : allProbes){
        if (aProbe.Area<(strongestArea*relativeAreaCutoff_))
          continue;
        if (aProbe.Area<(strongestArea*relativeFarAreaCutoff_) && (this.checkTimeSpace(strongestProbe,aProbe)>relativeFarAreaTimeSpace_))
          continue;
        Vector<CgProbe> finalProbes = new Vector<CgProbe>();
        finalProbes.add(aProbe);
        Hashtable<Integer,Vector<CgProbe>> isotopes = new Hashtable<Integer,Vector<CgProbe>>();
        isotopes.put(0, finalProbes);
        finalResults.put(hitNumber, isotopes);
        hitNumber++;
      }
    }catch (CgException | RulesException | NoRuleException | IOException | SpectrummillParserException | HydroxylationEncodingException | ChemicalFormulaException ex){
      ex.printStackTrace();
    }
    return finalResults;
  }
  
  

  /**
   * calculates a shotgun intensity
   * @param mz the m/z value of the analyte
   * @param mzTolerance the m/z tolerance for the extraction
   * @param charge the charge of the analyte
   * @param msLevel the MS-level
   * @return the detected VO
   * @throws CgException
   */
  private CgProbe calculateAShotgunIntensity(float mz, float mzTolerance, int charge, int msLevel) throws CgException{
    CgChromatogram cgChrom = readAChromatogram(mz, mzTolerance, mzTolerance, msLevel, 0f, 0);
    return calculateAShotgunIntensity(cgChrom, charge, msLevel, shotgunType_);
  }
  
  /**
   * calculates a shotgun intensity
   * @param cgChrom the chromatogram
   * @param charge the charge of the analyte
   * @param msLevel the MS-level
   * @param shotgunType the type of shotgun processing (i.e. SHOTGUN_TYPE_MEAN, SHOTGUN_TYPE_MEDIAN, or SHOTGUN_TYPE_SUM)
   * @return the detected VO
   */
  public static CgProbe calculateAShotgunIntensity(CgChromatogram cgChrom, int charge, int msLevel, int shotgunType) {
    CgProbe probe = new CgProbe(0,charge);
    probe.Area = calculateShotgunIntensity(cgChrom,shotgunType);
    probe.AreaStatus = CgAreaStatus.OK;
    probe.LowerMzBand = cgChrom.LowerMzBand;
    probe.UpperMzBand = cgChrom.UpperMzBand;
    probe.LowerValley = cgChrom.Value[0][0];
    probe.UpperValley = cgChrom.Value[cgChrom.Value.length-1][0];
    probe.Peak = (probe.LowerValley+probe.UpperValley)/2f;
    probe.Mz = cgChrom.Mz;
    return probe;
  }
  
  @SuppressWarnings("unchecked")
  private Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processByMzProbabsAndPossibleRetentionTime(float mz, int charge, float retentionTime, Vector<Float> retentionTimes, float prevTimeTolerance, float afterTimeTolerance, int timeType, Vector<Double>probabs, Vector<Double>possibleProbabs,int msLevel, boolean negative) throws CgException{
    coarseChromMzTolerance_ = LipidomicsConstants.getCoarseChromMzTolerance(mz);
    useSameCgHashFor3D_ = true;
    this.initCacheHashes();
    float massToAdd = LipidomicsConstants.getNeutronMass();
    if (negative) massToAdd*=-1;
//    this.sameCgHash_ = new Vector<CgProbe>();
//    this.same3DProbes_ = new Vector<Probe3D>();
//    this.same3DProbes_ = new Vector<Probe3D>();
    // this hash stores the chromatograms for every isotope
    chromHash_ = new Hashtable<Integer,LipidomicsChromatogram>();
    // this hash stores not allowed isotopic probes of a previous chromatogram
    prevIsoProbes_ = new Hashtable<Float,CgProbe>();

    // here a multiplication factor is calculated for the least intensity to declare
    // a peak from a previous chromatogram as valid
    double possibleIsoLowerMinimum = 2d/3d;
    if (possibleProbabs.size()>1){
      // exclude Areas which are probably from a different isotopic pattern
      possibleIsoLowerMinimum = 100;
      for (int i=0; i!=(possibleProbabs.size()-1);i++){
        double ratio = (possibleProbabs.get(i)/possibleProbabs.get(i+1));
        //this is if I search for negative isotopes
        //if (ratio<1) ratio = 1d/ratio;
        ratio *= (2d/3d);
        if (ratio<possibleIsoLowerMinimum)
          possibleIsoLowerMinimum = ratio;
      }
    }
    // this is the the -1 chromatogram to detect if the peak comes from a different isotope
    LipidomicsChromatogram previousIsoChrom = null;
    // this is the current chromatogram (isotope 0)
    LipidomicsChromatogram chrom = null;
    if (this.removeIfOtherIsotopePresent_){
      previousIsoChrom = new LipidomicsChromatogram(readAChromatogram(mz-massToAdd/(float)charge, coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
      previousIsoChrom.GetMaximumAndAverage();
      chrom = new LipidomicsChromatogram(readAChromatogram(mz, coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
      chrom.GetMaximumAndAverage();
    }else{
      // if we do not remove peaks from a previous isotopic chromatogram, we use a fake chromatogram with no intensity values
      CgChromatogram cgChrom = readAChromatogram(mz, coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_);
      chrom = new LipidomicsChromatogram(cgChrom);
      chrom.GetMaximumAndAverage();
      CgChromatogram prevChrom = new CgChromatogram(cgChrom.ScanCount);
      for (int i=0; i!=cgChrom.ScanCount; i++){
        prevChrom.Value[i][0] = cgChrom.Value[i][0];
        prevChrom.Value[i][1] = 0;
        prevChrom.Value[i][2] = 0;
        prevChrom.Value[i][3] = 0;
      }
      previousIsoChrom = new LipidomicsChromatogram(prevChrom);
      previousIsoChrom.GetMaximumAndAverage();
    }

    if (useNoiseCutoff_)
      this.calculateNoiseCutoffValue(chrom);
    if (minimumRelativeIntensity_!=null)
      this.calculateMinIntThreshold(chrom);
    chromHash_.put(0, chrom);
    Vector<Vector<CgProbe>> both = scanForHits(chrom,previousIsoChrom,possibleIsoLowerMinimum,charge,msLevel,retentionTime, retentionTimes, prevTimeTolerance, afterTimeTolerance, timeType);
    // stores all probes that originate from a different isotope
    Vector<CgProbe> otherIsoProbes1 = both.get(1);
    // stores all valid probes returned by the detection
    Vector<CgProbe> allProbes = both.get(0);
    // The first Vector still contains duplicates;
    // this method removes the duplicates
    Vector<CgProbe> probes1 = new Vector<CgProbe>();
    float biggestArea = 0;
    CgProbe highestProbe = null;
    for (int i=0;i!=allProbes.size();i++){   
      CgProbe aProbe = allProbes.get(i);
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        probes1.add(aProbe);
        if (aProbe.Area>biggestArea){
          biggestArea = aProbe.Area;
          highestProbe = aProbe;
        }  
      }  
    }
//    if (verbose_) System.out.println("probes1: "+probes1.size());
    // this removes very small peak, which are uninteresting for further processing
    Vector<CgProbe> probes2 = new Vector<CgProbe>();
    for (CgProbe aProbe : probes1){
      if (aProbe.Area>(biggestArea/peakDiscardingAreaFactor_))
          probes2.add(aProbe);
    }
//    if (verbose_) System.out.println("probes2: "+probes2.size());
    //checks if there is an iso probe between the two probes();
    Vector<CgProbe> probes3 = new Vector<CgProbe>();
    for (int i=0 ; i!=probes2.size();i++){
      CgProbe probe1 = probes2.get(i);
      boolean isFromOtherIso = false;
      for (int j=0;j!=probes2.size();j++){
        if (i!=j){
          CgProbe probe2 = probes2.get(j);
          // checks first if the peaks are close to one another and then if there is an isotopic probe in between them
          if (this.checkTimeSpace(probe1, probe2)<isotopeInBetweenTime_&&(this.otherIsoProbeInBetween(probe1,probe2,new Vector<CgProbe>(prevIsoProbes_.values())))){
            for (CgProbe otherIsoProbe:prevIsoProbes_.values()){
              boolean isoInBetween = false;
              // if there is an isotopic probe in between every one of them has to be checked which one is 
              // the right one; there can be several ones
              if (probe1.Peak>probe2.Peak){
                if (probe2.Peak<otherIsoProbe.Peak&&otherIsoProbe.Peak<probe1.Peak)
                  isoInBetween = true;
              }else{
                if (probe1.Peak<otherIsoProbe.Peak&&otherIsoProbe.Peak<probe2.Peak)
                  isoInBetween = true;      
              }
              if (isoInBetween){
                // if the area is 3x smaller take it is the one from the other isotope;
                if (probe1.Area<(probe2.Area)/isoInBetweenAreaFactor_){
                  isFromOtherIso = true;
                // if the area is comparable to the other one it has to be checked which one covers more parts of the isotope  
                }else if (probe1.Area<(probe2.Area)*isoInBetweenAreaFactor_){
                  float totalSpan = otherIsoProbe.UpperValley-otherIsoProbe.LowerValley;
                  if (probe1.Peak>probe2.Peak){
                    if (probe1.LowerValley<otherIsoProbe.UpperValley && probe2.UpperValley>otherIsoProbe.LowerValley){
                      float ratio1 = (otherIsoProbe.UpperValley-probe1.LowerValley)/totalSpan;
                      float ratio2 = (probe2.UpperValley-otherIsoProbe.LowerValley)/totalSpan;
//                    System.out.println(ratio1+";"+ratio2);
                      // if this peak covers more percent of the other isotope it is regarded as the one from the other isotope
                      if (ratio1>ratio2){
                        isFromOtherIso = true;
                      }
                    }
                  }else{
                    if (probe2.LowerValley<otherIsoProbe.UpperValley && probe1.UpperValley>otherIsoProbe.LowerValley){
                      float ratio1 = (probe1.UpperValley-otherIsoProbe.LowerValley)/totalSpan;
                      float ratio2 = (otherIsoProbe.UpperValley-probe2.LowerValley)/totalSpan;                  
//                    System.out.println(ratio1+";"+ratio2);
                      // if this peak covers more percent of the other isotope it is regarded as the one from the other isotope
                      if (ratio1>ratio2){
                        isFromOtherIso = true;
                      }
                  }                  
                }
                }
              }
            }
          }
        }
      }
//      System.out.println("isFromOtherIso: "+isFromOtherIso+";"+probe1.Area);
      // if the peak belongs to another isotope it is flagged and added to the otherIsoProbes Vector if not it comes to probes3
      if (isFromOtherIso){
        probe1.AreaStatus = CgAreaStatus.OtherIso;
        otherIsoProbes1.add(probe1);
      }else{
        probes3.add(probe1);
      }
    }
//    if (verbose_) System.out.println("probes3: "+probes3.size());

//    System.out.println("otherIsoProbes1: "+otherIsoProbes1.size());
    // isotopic probes that are far away from the accepted probes are out of interest
    // hence they are filtered out and just the ones near accepted probes are kept
    Vector<CgProbe> probes4 = this.uniteTwinPeaks(probes3,chrom,msLevel);
    Vector<CgProbe> otherIsoProbes2 = new Vector<CgProbe>();
    for (CgProbe otherIsoProbe: otherIsoProbes1){
      if (otherIsoProbe.AreaStatus == CgAreaStatus.OtherIso && this.checkNearNormalProbe(otherIsoProbe,probes4,null)){
        otherIsoProbes2.add(otherIsoProbe);
      }  
    }
//    if (verbose_) System.out.println("probes4: "+probes4.size());

//    System.out.println("otherIsoProbes2: "+otherIsoProbes2.size());
    
    boolean[] isOfInterest = new boolean[probes4.size()];
//    double mainPeak = 1;
    Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicAreas = new Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>>();
    for (int j=0; j!=probes4.size();j++){
      isOfInterest[j] = true;
      Vector<CgProbe> areas = new Vector<CgProbe>();
      areas.add(probes4.get(j));
      Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneIso = new Hashtable<Integer,Vector<CgProbe>>();
      isoAreasOfOneIso.put(0, areas);
      isotopicAreas.put(j, isoAreasOfOneIso);
    }
//    if (probabs.size()>0)
//      mainPeak = probabs.get(0);
    Vector<CgProbe> currentIsoProbes = new Vector<CgProbe>(otherIsoProbes2);
    Vector<CgProbe> recentIsoProbes;
    // now the other isotopic peaks of one series that have to be fulfilled are calculated
    // if one probe does not fulfil the theoretical isotopic area-ratio, it is discarded (isOfInterest[j]=false)
    // it iterates first over the isotopic series starting with isotope 1
    for (int i=0;i!=probabs.size();i++){
//      previousIsoChrom = chrom;
      if (i!=0){
        recentIsoProbes = currentIsoProbes;
        currentIsoProbes = new Vector<CgProbe>();
        // the chromatogram of the specific series is extracted
        chrom = new LipidomicsChromatogram(readAChromatogram(mz+i*massToAdd/(float)charge, coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        chrom.GetMaximumAndAverage();
        chromHash_.put(i, chrom);
        // the multiplicationFactor defines theoretical ratio between the 0-isotope and the one try to calculate
        double multiplicationFactor = probabs.get(i)/probabs.get(0);
        // now calculate the isotopic area for each probe that is still OK and discard it if it does not fulfil the ratio
        for (int j=0; j!=probes4.size();j++){
          if (isOfInterest[j]){
            Probe3D aProbe = (Probe3D)probes4.get(j);
            Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneProbe = new Hashtable<Integer,Vector<CgProbe>>();
            if (isotopicAreas.containsKey(j))
              isoAreasOfOneProbe = isotopicAreas.get(j);
            // this method takes care if all of the criteria are fulfilled and returns the corresponding isotopic area
            @SuppressWarnings("rawtypes")
            Vector result = calculateIsoProbesWithinZeroIsotopeBoundaries(chrom,aProbe,probes4,multiplicationFactor,recentIsoProbes,probabs.get(1)/probabs.get(0),false,charge,msLevel);
            CgProbe isoProbe = (CgProbe)result.get(0);
            currentIsoProbes.addAll((Vector<CgProbe>)result.get(1));
//            if (isoProbe!=null) System.out.println(isoProbe.AreaStatus);
            // if there is a returned probe and it is OK than it is valid, otherwise discard it (isOfInterest[j]=false)
            if (isoProbe!=null&&isoProbe.AreaStatus==CgAreaStatus.OK){
              isoProbe.isotopeNumber = i;
              Vector<CgProbe> isoAreas = new Vector<CgProbe>();
              isoAreas.add(isoProbe);
//              System.out.println(i+": "+isoProbe.Area/aProbe.Area+";"+multiplicationFactor+";"+(isoProbe.Area/aProbe.Area)/multiplicationFactor);
              isoAreasOfOneProbe.put(i,isoAreas);
              isotopicAreas.put(j, isoAreasOfOneProbe);
            }else{
              isOfInterest[j] = false;
              isotopicAreas.remove(j);
            }            
          }
        }
        currentIsoProbes = this.curateCurrentIsoProbes(currentIsoProbes,i,isotopicAreas,isOfInterest);
      }
    }
//    int countOK = 0;
//    for (int i=0;i!=isOfInterest.length;i++){
//      if (isOfInterest[i]){
//        countOK++;
//      }else{
//
//      }
//    }
//    if (verbose_)System.out.println("CountOK: "+countOK);
    // remove the small Peaks
    biggestArea = 0;
    highestProbe = null;
    for (int i=0;i!=probes4.size();i++){
      if (isOfInterest[i]){
        CgProbe aProbe = probes4.get(i);
        if (aProbe.AreaStatus==CgAreaStatus.OK){
          if (aProbe.Area>biggestArea){
            biggestArea = aProbe.Area;
            highestProbe = aProbe;
          }
        }
      }
    }
    // isotopic probes that are far away from the accepted probes are out of interest
    // hence they are filtered out and just the ones near accepted probes are kept
    Vector<CgProbe> otherIsoProbes3 = new Vector<CgProbe>();
    for (CgProbe otherIsoProbe: otherIsoProbes2){
      if (otherIsoProbe.AreaStatus == CgAreaStatus.OtherIso && this.checkNearNormalProbe(otherIsoProbe,probes4,isOfInterest)){
        otherIsoProbes3.add(otherIsoProbe);
      }  
    }
//    System.out.println("otherIsoProbes3: "+otherIsoProbes3.size());
    // here another filter is applied to remove small peaks
    // first, the found area has to have at least 1% of the highest found are (the ones from a different isotope are not counted)
    // second, if the area is farer away from the highest probe it has to have at least 10% of its area; the reason for this
    // is that very close ones could be part of twin-peaks and should not be removed
    // third, if there is an other-isotope-probe between the highest probe and this probe the area has to be at least half of 
    // the highest area, because then we cannot really distinguish anymore which one belongs to the other isotope
    for (int i=0;i!=probes4.size();i++){
      if (isOfInterest[i]){
        CgProbe aProbe = probes4.get(i);
        if (aProbe.AreaStatus==CgAreaStatus.OK){
          if (aProbe.Area>(biggestArea*relativeAreaCutoff_)){
//            System.out.println("NormalProbeArea: "+aProbe.Area);
            if (aProbe.Area<(biggestArea*relativeFarAreaCutoff_) && (this.checkTimeSpace(highestProbe,aProbe)>relativeFarAreaTimeSpace_)){
              isOfInterest[i] = false;
              isotopicAreas.remove(i);
            }
            if (aProbe.Area<(biggestArea*relativeIsoInBetweenCutoff_) && (otherIsoProbeInBetween(highestProbe,aProbe,otherIsoProbes3))){
              isOfInterest[i] = false;
              isotopicAreas.remove(i);
            }
          }else{
            isOfInterest[i] = false;
            isotopicAreas.remove(i);
          }
        }
      }
    }
//    countOK = 0;
//    for (int i=0;i!=isOfInterest.length;i++){
//      if (isOfInterest[i])
//        countOK++;
//    }
//    System.out.println("CountOK: "+countOK);
    // here is the same procedure as for the "must" isotopic probes 
    // the only difference: if there is no isotope has been found the probe is not removed
    if (possibleProbabs.size()>probabs.size()){
//      mainPeak = possibleProbabs.get(0);
      if (probabs.size()==0)
        recentIsoProbes = new Vector<CgProbe>(otherIsoProbes3);
      for (int i=probabs.size(); i!=possibleProbabs.size();i++){
        if (i!=0){
          recentIsoProbes = currentIsoProbes;
          currentIsoProbes = new Vector<CgProbe>();
          chrom = new LipidomicsChromatogram(readAChromatogram(mz+(i*massToAdd)/(float)charge, coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
          chrom.GetMaximumAndAverage();
          chromHash_.put(i, chrom);
          double multiplicationFactor = possibleProbabs.get(i)/possibleProbabs.get(0);
          for (int j=0; j!=probes4.size();j++){
            if (isOfInterest[j]){
              Probe3D aProbe = (Probe3D)probes4.get(j);
              Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneProbe = new Hashtable<Integer,Vector<CgProbe>>();
              if (isotopicAreas.containsKey(j))
                isoAreasOfOneProbe = isotopicAreas.get(j);
              @SuppressWarnings("rawtypes")
              Vector result = calculateIsoProbesWithinZeroIsotopeBoundaries(chrom,aProbe,probes4,multiplicationFactor,recentIsoProbes,possibleProbabs.get(1)/possibleProbabs.get(0),false,charge,msLevel);
              CgProbe isoProbe = (CgProbe)result.get(0);
              currentIsoProbes.addAll((Vector<CgProbe>)result.get(1));
              if (isoProbe!=null&&isoProbe.AreaStatus==CgAreaStatus.OK){
                isoProbe.isotopeNumber = i;
                Vector<CgProbe> isoAreas = new Vector<CgProbe>();
//                System.out.println(i+": "+isoProbe.Area/aProbe.Area+";"+multiplicationFactor+";"+(isoProbe.Area/aProbe.Area)/multiplicationFactor);
                isoAreas.add(isoProbe);
                isoAreasOfOneProbe.put(i,isoAreas);
                isotopicAreas.put(j, isoAreasOfOneProbe);
              }       
            }
          }
          currentIsoProbes = this.curateCurrentIsoProbes(currentIsoProbes,i,isotopicAreas,isOfInterest);
        }
      }
    }
    // now the hash-table is a little bit changed; this hash contains first the isotope number
    // and then a vector of the added probes; the old one was sorted based on the probe, then isotope and then values
    isotopicAreas = uniteClosePeaks(uniteClosePeaks(isotopicAreas,chromHash_,msLevel),chromHash_,msLevel);
    this.releaseCacheHashes();
    return isotopicAreas;
  }
  
  private Vector<Vector<CgProbe>> scanForHits(LipidomicsChromatogram chrom, LipidomicsChromatogram previousIsoChrom, double possibleIsoLowerMinimum, int charge, int msLevel, float retentionTime, Vector<Float> retentionTimes, float prevTimeTolerance, float afterTimeTolerance, int timeType) throws CgException{
    Vector<Vector<CgProbe>> both = new Vector<Vector<CgProbe>>();
    Vector<CgProbe> allProbes = new Vector<CgProbe>();
    Vector<CgProbe> otherIsoProbes1 = new Vector<CgProbe>();
    // the current probe of interest
    Probe3D probe = null;
    int currentScan = 0;
    Probe3D currentValidProbe = null;
    if (retentionTimes!=null && retentionTimes.size()>0){
      float prevTime = prevTimeTolerance;
      float afterTime = afterTimeTolerance;
      if (timeType==LipidomicsDefines.MINUTES){
        prevTime = prevTime*60f;
        afterTime = afterTime*60f;
      }
      //first check if there are peaks at the found spectra
      for (float retTime : retentionTimes){
        boolean isCoveredByProbe = false;
        for (CgProbe aProbe : allProbes){
          if (aProbe.LowerValley<retTime && retTime<aProbe.UpperValley){
            isCoveredByProbe = true;
            break;
          }
        }
        if (isCoveredByProbe) continue;
        for (CgProbe aProbe : otherIsoProbes1){
          if (aProbe.LowerValley<retTime && retTime<aProbe.UpperValley){
            isCoveredByProbe = true;
            break;
          }
        }
        if (isCoveredByProbe) continue;
        int mainScan = LipidomicsAnalyzer.findIndexByTime(retTime, chrom);
        probe = this.calculateAreaWithoutPrevIso(chrom,previousIsoChrom, mainScan,possibleIsoLowerMinimum,true,charge,msLevel);
        if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
          // checks if such a peak is already there
          this.checkDuplicate(allProbes,probe);
          if (probe.AreaStatus==CgAreaStatus.OK)
            allProbes.add(probe);
          // if the peak belongs to another isotope, the peak is put into a separate vector    
        }else if (probe!=null&&probe.AreaStatus==CgAreaStatus.OtherIso){
          this.checkDuplicate(otherIsoProbes1,probe);
          if (probe.AreaStatus==CgAreaStatus.OtherIso)
            otherIsoProbes1.add(probe);
        }        
      }
      //second check if there are regions left that where not covered by peaks
      Vector<Range> rangesToCover = new Vector<Range>();
      for (float retTime : retentionTimes){
        Range range = new Range(retTime-prevTime,retTime+afterTime);
        boolean overlap = false;
        for (int i=0; i!=rangesToCover.size(); i++){
          Range exRange = rangesToCover.get(i);
          if (range.overlap(exRange)){
            exRange.combine(range);
            overlap = true;
            break;
          }
        }
        if (!overlap) rangesToCover.add(range);
      }
      for (CgProbe aProbe : allProbes){
        Vector<Range> reducedCover = new Vector<Range>();
        for (Range range : rangesToCover){
          // range is covered by found probe -> remove it
          reducedCover.addAll(Range.reduce(range, new Range(aProbe.LowerValley,aProbe.UpperValley)));          
        }
        rangesToCover = reducedCover;
      }
      for (CgProbe aProbe : otherIsoProbes1){
        Vector<Range> reducedCover = new Vector<Range>();
        for (Range range : rangesToCover){
          // range is covered by found probe -> remove it
          reducedCover.addAll(Range.reduce(range, new Range(aProbe.LowerValley,aProbe.UpperValley)));          
        }
        rangesToCover = reducedCover;
      }
      // now scan the ranges
      for (Range range : rangesToCover){
        currentScan = LipidomicsAnalyzer.findIndexByTime(range.getStart(), chrom);
        // the highest scan when looking in the time range after the main scan
        int stopScan = LipidomicsAnalyzer.findIndexByTime(range.getStop(), chrom)+1;
        while (currentScan<chrom.ScanCount&&currentScan<stopScan){
          // calculates a peak/probe for the scan and checks if it is from a different isotope 
          probe = this.calculateAreaWithoutPrevIso(chrom,previousIsoChrom, currentScan,possibleIsoLowerMinimum,true,charge,msLevel);
          // if the probe is OK it is added to all probes and the upper scan index
          // is set to the upper border of the peak; this reduces the amount of unnecessary calculations
          if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
            this.checkDuplicate(allProbes,probe);
            if (probe.AreaStatus==CgAreaStatus.OK)
              allProbes.add(probe);
            if (currentValidProbe!=null){
              int upScan = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.UpperValley,chrom);
              if (upScan>currentScan)
                currentScan = upScan;
            }
            // if the peak belongs to another isotope, the peak is put into a separate vector  
          }else if (probe!=null && probe.AreaStatus==CgAreaStatus.OtherIso){
            this.checkDuplicate(otherIsoProbes1,probe);
            if (probe.AreaStatus==CgAreaStatus.OtherIso)
              otherIsoProbes1.add(probe);
          }
          currentScan +=scanStep_;
          //this is for duplicate removal
          Vector<CgProbe> probes = new Vector<CgProbe>(allProbes);
          allProbes = new Vector<CgProbe>();
          for (CgProbe another : probes){
            if (another.AreaStatus != CgAreaStatus.Duplicate) allProbes.add(another);
          }
        }
      }
      
    } else {
      int stopScan = Integer.MAX_VALUE;

      // if just a specific retention time space is scanned we start at a specific
      // scan go down some scans (settable amount of minutes);
      // for the positive time range we use the same method like we do not specify 
      // a retention time; we just set the variable currentScan to the highest scan number
      // so far and start there; the no retention time scan simply starts at zero and 
      // proceeds then one by one step
      if (retentionTime>0 && prevTimeTolerance>0 && afterTimeTolerance>0){
        float mainRetTime = retentionTime;
        float startRetTime = mainRetTime-prevTimeTolerance;
        if (startRetTime<0)
          startRetTime = 0f;
        float stopRetTime = mainRetTime+afterTimeTolerance;
        if (timeType == LipidomicsDefines.MINUTES){
          mainRetTime = mainRetTime*60f;
          startRetTime = startRetTime*60f;
          stopRetTime = stopRetTime*60f;
        }
        // the scan where to start the search for the peak
        int mainScan = LipidomicsAnalyzer.findIndexByTime(mainRetTime, chrom);
        // the lowest scan when looking in the time range before the main scan
        int startScan = LipidomicsAnalyzer.findIndexByTime(startRetTime, chrom);
        // the highest scan when looking in the time range after the main scan
        stopScan = LipidomicsAnalyzer.findIndexByTime(stopRetTime, chrom)+1;
        // calculates a peak/probe for the scan and checks if it is from a different isotope 
        probe = this.calculateAreaWithoutPrevIso(chrom,previousIsoChrom, mainScan,possibleIsoLowerMinimum,true,charge,msLevel);
        currentScan = mainScan;
        // if the probe is OK it is added to all probes and the lower scan index
        // is set to the lower border of the peak; this reduces the amount of unnecessary calculations
        if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
          allProbes.add(probe);
          currentValidProbe = probe;
          int lowScan = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.LowerValley,chrom);
          if (lowScan<currentScan)
            currentScan = lowScan;
          // if the peak belongs to another isotope, the peak is put into a separate vector  
        }else if (probe!=null&&probe.AreaStatus==CgAreaStatus.OtherIso){
          otherIsoProbes1.add(probe);
        }
        currentScan--;
        // now we go to the left side of the peak and look within the specified range for other possible peaks
        while (currentScan>=startScan){
          // calculates a peak/probe for the scan and checks if it is from a different isotope 
          probe = this.calculateAreaWithoutPrevIso(chrom,previousIsoChrom, currentScan,possibleIsoLowerMinimum,true,charge,msLevel);
          // if the probe is OK it is added to all probes and the lower scan index
          // is set to the lower border of the peak; this reduces the amount of unnecessary calculations
          if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
            // checks if such a peak is already there
            this.checkDuplicate(allProbes,probe);
            if (probe.AreaStatus==CgAreaStatus.OK)
              allProbes.add(probe);
            if (currentValidProbe!=null){
              int lowScan = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.LowerValley,chrom);
              if (lowScan<currentScan)
                currentScan = lowScan;
            }
            // if the peak belongs to another isotope, the peak is put into a separate vector    
          }else if (probe!=null&&probe.AreaStatus==CgAreaStatus.OtherIso){
            this.checkDuplicate(otherIsoProbes1,probe);
            if (probe.AreaStatus==CgAreaStatus.OtherIso)
              otherIsoProbes1.add(probe);
          }
          Vector<CgProbe> probes = new Vector<CgProbe>(allProbes);
          allProbes = new Vector<CgProbe>();
          for (CgProbe another : probes){
            if (another.AreaStatus != CgAreaStatus.Duplicate) allProbes.add(another);
          }
          currentScan -= scanStep_;
        }
        currentScan = mainScan;
        // this checks if the current scan for the right side direction 
        // should be put at the top border of the peak or if the start scan should be used 
        if (allProbes.size()>0){
          CgProbe highestTimeProbe = allProbes.get(0);
          if (currentValidProbe!=null){
            int highesTimeScan = LipidomicsAnalyzer.findIndexByTime(highestTimeProbe.UpperValley,chrom);
            if (highesTimeScan>currentScan)
              currentScan = highesTimeScan;
          }
        }
        currentScan ++;
      }
      // this procedure goes from a start scan (can be 0) to the right until it reaches
      // the end or a defined border
      while (currentScan<chrom.ScanCount&&currentScan<stopScan){
        // calculates a peak/probe for the scan and checks if it is from a different isotope 
        probe = this.calculateAreaWithoutPrevIso(chrom,previousIsoChrom, currentScan,possibleIsoLowerMinimum,true,charge,msLevel);
        // if the probe is OK it is added to all probes and the upper scan index
        // is set to the upper border of the peak; this reduces the amount of unnecessary calculations
        if (probe!=null&&probe.AreaStatus==CgAreaStatus.OK){
          this.checkDuplicate(allProbes,probe);
          if (probe.AreaStatus==CgAreaStatus.OK)
            allProbes.add(probe);
          if (currentValidProbe!=null){
            int upScan = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.UpperValley,chrom);
            if (upScan>currentScan)
              currentScan = upScan;
          }
          // if the peak belongs to another isotope, the peak is put into a separate vector  
        }else if (probe!=null && probe.AreaStatus==CgAreaStatus.OtherIso){
          this.checkDuplicate(otherIsoProbes1,probe);
          if (probe.AreaStatus==CgAreaStatus.OtherIso)
            otherIsoProbes1.add(probe);
        }
        currentScan +=scanStep_;
        //this is for duplicate removal
        Vector<CgProbe> probes = new Vector<CgProbe>(allProbes);
        allProbes = new Vector<CgProbe>();
        for (CgProbe another : probes){
          if (another.AreaStatus != CgAreaStatus.Duplicate) allProbes.add(another);
        }
      }
    }
    both.add(allProbes);
    both.add(otherIsoProbes1);
    return both;
  }

  @SuppressWarnings("unchecked")
  public Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processByMzAndRetentionTime(float mz, int charge, //float mzTolerance, 
      float retentionTime, float prevTimeTolerance, float afterTimeTolerance, int timeType, Vector<Double>probabs,
      Vector<Double>possibleProbabs, int msLevel, boolean negative) throws CgException{
    this.coarseChromMzTolerance_ = LipidomicsConstants.getCoarseChromMzTolerance(mz);
    float massToAdd = LipidomicsConstants.getNeutronMass();
    if (negative) massToAdd*=-1;
//    return this.processByMzProbabsAndPossibleRetentionTime(mz, //mzTolerance, 
//        retentionTime, prevTimeTolerance, afterTimeTolerance, timeType, probabs, probabs);
    
    //m_params = new CgParameterSet(mz,molecule,mz,-1,-1,-1,-1,-1,-1,-1,-1);
    
//    Vector<Vector<CgProbe>> isotopicProbes = new Vector<Vector<CgProbe>>();
    Vector<CgProbe> allProbes = new Vector<CgProbe>();
    float mainRetTime = retentionTime;
    float startRetTime = mainRetTime-prevTimeTolerance;
    if (startRetTime<0)
      startRetTime = 0f;
    float stopRetTime = mainRetTime+afterTimeTolerance;
    if (timeType == LipidomicsDefines.MINUTES){
      mainRetTime = mainRetTime*60f;
      startRetTime = startRetTime*60f;
      stopRetTime = stopRetTime*60f;
    }
    LipidomicsChromatogram previousIsoChrom = new LipidomicsChromatogram(readAChromatogram(mz-massToAdd/(float)charge, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
    previousIsoChrom.GetMaximumAndAverage();
    double multiplicationFactor = 1;
    if (probabs.size()>1){
      // exclude Areas which are probably from a different isotopic pattern
      multiplicationFactor = (probabs.get(0)/probabs.get(1))*(2f/3f);
    }  

    LipidomicsChromatogram chrom = new LipidomicsChromatogram(readAChromatogram(mz, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
    int mainScan = LipidomicsAnalyzer.findIndexByTime(mainRetTime, chrom);
    int startScan = LipidomicsAnalyzer.findIndexByTime(startRetTime, chrom);
    int stopScan = LipidomicsAnalyzer.findIndexByTime(stopRetTime, chrom);
    chrom.GetMaximumAndAverage();
    @SuppressWarnings("rawtypes")
    Vector result = this.calculateOneAreaWithoutPrevIso(chrom,previousIsoChrom, mainScan,multiplicationFactor,new Vector<CgProbe>(),true,charge);
    Vector<CgProbe> theProbes = (Vector<CgProbe>)result.get(0);
    int lowestScan = (Integer)result.get(1);
    int highestScan = (Integer)result.get(2);
//    System.out.println("theProbes.size(): "+theProbes.size());
//    for (CgProbe probe: theProbes){
//      System.out.println(probe.Area+";"+probe.Peak+";"+probe.LowerValley+";"+probe.UpperValley+";"+probe.AreaStatus);
//    }

    
    
    //Since the returned probes start with the one with the lowest retention time, I have the start at the top of the list,
    //because I want to go downwards
    CgProbe currentValidProbe = null;
    for (int i=(theProbes.size()-1);i!=-1;i--){
      CgProbe anProbe = theProbes.get(i);
//      System.out.println("anProbe.Area: "+anProbe.Area);
      if (anProbe!=null&&anProbe.AreaStatus==CgAreaStatus.OK){
        allProbes.add(anProbe);
        currentValidProbe = anProbe;
      }
    }
    int currentScan = mainScan-1;
    if (currentValidProbe!=null){
      currentScan = lowestScan;
//      currentScan = this.findIndexByTime(currentValidProbe.LowerValley,chrom);
    }
    while (currentScan>=startScan){
      result = this.calculateOneAreaWithoutPrevIso(chrom,previousIsoChrom, currentScan,multiplicationFactor,new Vector<CgProbe>(),true,charge);
      theProbes = (Vector<CgProbe>)result.get(0);
      lowestScan = (Integer)result.get(1);
//      probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.EnhancedValleyMethod);
      //Since the returned probes start with the one with the lowest retention time, I have the start at the top of the list,
      //because I want to go downwards
//      System.out.println(currentScan+": "+theProbes.size());
      for (int i=(theProbes.size()-1);i!=-1;i--){
        CgProbe probe = theProbes.get(i);
        if (probe.AreaStatus==CgAreaStatus.OK){
          this.checkDuplicate(allProbes,probe);
          if (!probe.greedyProbe&&probe.AreaStatus==CgAreaStatus.Duplicate){
            probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.StandardValleyMethod,charge);
          }
          int isPossibleIsotope = 0;
          if (probe.greedyProbe)
            isPossibleIsotope = isFromOtherIsotope(probe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>());
          else
            isPossibleIsotope = isFromOtherIsotope(probe, LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>());
          if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP&&(!probe.greedyProbe)){
            CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
            greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(probe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,new Vector<CgProbe>(),charge);
  //          if (!isGreedyFromOtherIsotope(greedyProbe,probe,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>())){
//            if (isFromOtherIsotope(greedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom, multiplicationFactor,new Vector<CgProbe>())!= LipidomicsAnalyzer.ISOTOPE_PRESENT){
            if (greedyProbe!=null){
              greedyProbe.greedyProbe = true;
              probe = greedyProbe;
            }//else{
             // isPossibleIsotope=LipidomicsAnalyzer.ISOTOPE_PRESENT;
            //}
          }
          if (isPossibleIsotope!=LipidomicsAnalyzer.ISOTOPE_PRESENT){
              this.checkDuplicate(allProbes,probe);
            if (currentValidProbe==null||probe.AreaStatus==CgAreaStatus.OK){
              currentValidProbe = probe;
            }
            allProbes.add(probe);
          }  
//        int lowerValleyIndex = this.findIndexByTime(probe.LowerValley,chrom);
//        if (lowerValleyIndex<currentScan)
//          currentScan = lowerValleyIndex;
          if (lowestScan<currentScan)
            currentScan = lowestScan;
          
        }
      }
      currentScan--;
    }
    currentScan = mainScan+1;
    if (allProbes.size()>0){
//      int highesTimeScan = this.findIndexByTime(highestTimeProbe.UpperValley,chrom);
      int highesTimeScan = highestScan;
      if (highesTimeScan>currentScan)
        currentScan = highesTimeScan;
    }
    while (currentScan<=stopScan){
//      probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.EnhancedValleyMethod);
      result = this.calculateOneAreaWithoutPrevIso(chrom,previousIsoChrom, currentScan,multiplicationFactor,new Vector<CgProbe>(),true,charge);
      theProbes = (Vector<CgProbe>)result.get(0);
      highestScan = (Integer)result.get(2);
      for (int i=0;i!=theProbes.size();i++){
        CgProbe probe = theProbes.get(i);
        if (probe.AreaStatus==CgAreaStatus.OK){
          this.checkDuplicate(allProbes,probe);
          this.checkDuplicate(allProbes,probe);
          if (!probe.greedyProbe&&probe.AreaStatus==CgAreaStatus.Duplicate){
            probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.StandardValleyMethod,charge);
          }
          int isPossibleIsotope = 0;
          if (probe.greedyProbe)
            isPossibleIsotope = isFromOtherIsotope(probe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>());
          else
            isPossibleIsotope = isFromOtherIsotope(probe, LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>());//            if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP){
//              probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.GreedySteepnessReductionMethod);
//              probe.greedyProbe = true;
//            }
            if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP&&(!probe.greedyProbe)){
              CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
              greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(probe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,new Vector<CgProbe>(),charge);
//              if (!isGreedyFromOtherIsotope(greedyProbe,probe,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>())){
//              if (isFromOtherIsotope(greedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom, multiplicationFactor,new Vector<CgProbe>())!= LipidomicsAnalyzer.ISOTOPE_PRESENT){
              if (greedyProbe!=null){
                greedyProbe.greedyProbe = true;
                probe = greedyProbe;    
              }//else{
               // isPossibleIsotope=LipidomicsAnalyzer.ISOTOPE_PRESENT;
              //}
            }
            if (isPossibleIsotope!=LipidomicsAnalyzer.ISOTOPE_PRESENT){

              this.checkDuplicate(allProbes,probe);
              if (currentValidProbe==null||probe.AreaStatus==CgAreaStatus.OK){
                currentValidProbe = probe;
              }
              allProbes.add(probe);
            }
            if (highestScan>currentScan)
              currentScan = highestScan;
//            if (upperValleyIndex>currentScan)
//              currentScan = upperValleyIndex;
        }
      }
      currentScan++;
    }
//    System.out.println("allProbes.size(): "+allProbes.size());
    Vector<CgProbe> probes1 = new Vector<CgProbe>();
    float biggestArea = 0;
    CgProbe highestProbe = null;
    for (int i=0;i!=allProbes.size();i++){
      CgProbe aProbe = allProbes.get(i);
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        probes1.add(aProbe);
        if (aProbe.Area>biggestArea){
          biggestArea = aProbe.Area;
          highestProbe = aProbe;
        }  
      }  
    }
//    Vector<CgProbe> probes2 = new Vector<CgProbe>();
//    for (CgProbe aProbe : probes1){
//      if (aProbe.Area>(biggestArea/200)){
//        if (aProbe.Area>(biggestArea/10)&&this.checkTimeSpace(highestProbe,aProbe)<30)
//          probes2.add(aProbe);
//      }  
//    }
//    Vector<CgProbe> probes = probes2;
    Vector<CgProbe> probes = probes1;
//    if (probabs.size()>1){
//
//      // calculate the areas for the isotopic pattern
      double mainPeak = 1;
      boolean[] isOfInterest = new boolean[probes.size()];
//      if (probes.size()>0)
//        isotopicProbes.add(probes);
      for (int i=0; i!=probes.size();i++){
        isOfInterest[i] = true;
      }
//      for (int i=0;i!=probabs.size();i++){
//        previousIsoChrom = chrom;
//        if (i==0){
//          mainPeak = probabs.get(i);
//        }else{
//          double relativeValue = probabs.get(i)/mainPeak;
//          chrom = this.readAChromatogram(mz+i*NEUTRON_MASS, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, chromSmoothRange_,chromSmoothRepeats_);
//          chrom.GetMaximumAndAverage();
//          Vector<CgProbe> oneisotopeProbes = new Vector<CgProbe>();
//          for (int j=0; j!=probes.size();j++){
//            if (isOfInterest[j]){
//              CgProbe aProbe = probes.get(j);
//              int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, chrom);
//              result = this.calculateOneAreaWithoutPrevIso(chrom,previousIsoChrom, scan,multiplicationFactor,isotopicProbes.get(isotopicProbes.size()-1),false);
//              double theoreticalAreaValue =  aProbe.Area*relativeValue;
//              theProbes = (Vector<CgProbe>)result.get(0);
//              boolean isOneOfInterest = false;
//              for (CgProbe isotopeProbe : theProbes){
//                if (isotopeProbe.AreaStatus==CgAreaStatus.OK){
//                  CgProbe isotopeProbeStandard = LipidomicsAnalyzer.calculateOneArea(chrom, scan, LipidomicsDefines.StandardValleyMethod) ;
//                  CgProbe isotopeProbeGreedy = LipidomicsAnalyzer.calculateOneArea(chrom, scan, LipidomicsDefines.GreedySteepnessReductionMethod);
//                  isotopeProbeGreedy.greedyProbe = true;
//                  double[] isotopeRegionAreas = this.calculateLowerUpperArea(i, theoreticalAreaValue, aProbe.greedyProbe);
//                  double lowerArea = isotopeRegionAreas[0];
//                  double upperArea = isotopeRegionAreas[1];
//                  if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea) || (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
//                      (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
//                    CgProbe rightProbe = isotopeProbe;
//                    if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
//                        (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)||
//                        (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
//                      if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)&&
//                          isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
//                        rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard,isotopeProbeGreedy);
//                      }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
//                        rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard);
//                      }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
//                        rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeGreedy);
//                      } else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
//                        rightProbe = this.checkCloserToExpected(aProbe,isotopeProbeStandard,isotopeProbeGreedy);
//                      }
//                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea){
//                      rightProbe = isotopeProbe;
//                    }else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
//                      rightProbe = isotopeProbeStandard;
//                    }else if (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
//                      rightProbe = isotopeProbeGreedy;
//                    }
//                    this.checkDuplicate(oneisotopeProbes, rightProbe);
//                    oneisotopeProbes.add(rightProbe);
//
//                    isOneOfInterest = true;
//                  }
//                }
//              }
//              if (!isOneOfInterest)
//                isOfInterest[j] = false;
//            }
//          }
//          if (oneisotopeProbes.size()>0){
//            Vector<CgProbe> oneisotopeProbeOK = new Vector<CgProbe>();
//            for (CgProbe aProbe: oneisotopeProbes){
//              if (aProbe.AreaStatus==CgAreaStatus.OK){
//                oneisotopeProbeOK.add(aProbe);
//              }
//            }
//            isotopicProbes.add(oneisotopeProbeOK);
//          }
//        }
//      }
//      
//      
//    }else{
//      probes = probes2;
//      if (probes.size()>0)
//        isotopicProbes.add(probes);
//    }
//    return isotopicProbes;
    Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicAreas = new Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>>();
    for (int j=0; j!=probes.size();j++){
      if (isOfInterest[j]){
        Vector<CgProbe> areas = new Vector<CgProbe>();
        areas.add(probes.get(j));
        Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneIso = new Hashtable<Integer,Vector<CgProbe>>();
        isoAreasOfOneIso.put(0, areas);
        isotopicAreas.put(j, isoAreasOfOneIso);
      }
    }
    for (int i=0;i!=probabs.size();i++){
      previousIsoChrom = chrom;
      if (i==0){
        mainPeak = probabs.get(i);
      }else{
        double relativeValue = probabs.get(i)/mainPeak;
        chrom = new LipidomicsChromatogram(readAChromatogram(mz+(i*massToAdd)/(float)charge, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        chrom.GetMaximumAndAverage();
        Vector<CgProbe> allowedPrevIsos = this.getAllowedPrevIsos(isotopicAreas,(i-1),isOfInterest);
        for (int j=0; j!=probes.size();j++){
          if (isOfInterest[j]){
            CgProbe aProbe = probes.get(j);
            Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneProbe = new Hashtable<Integer,Vector<CgProbe>>();
            if (isotopicAreas.containsKey(j))
              isoAreasOfOneProbe = isotopicAreas.get(j);
            int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, chrom);
            result = this.calculateOneAreaWithoutPrevIso(chrom,previousIsoChrom, scan,multiplicationFactor,allowedPrevIsos,false,charge);
            theProbes = (Vector<CgProbe>)result.get(0);
            boolean isOneOfInterest = false;
            double theoreticalAreaValue =  aProbe.Area*relativeValue;
            Vector<CgProbe> oneisotopeProbes = new Vector<CgProbe>();
            for (CgProbe isotopeProbe : theProbes){
              isotopeProbe.isotopeNumber = i;
              if (isotopeProbe.AreaStatus==CgAreaStatus.OK){
                CgProbe isotopeProbeStandard = LipidomicsAnalyzer.calculateOneArea(chrom, scan, LipidomicsDefines.StandardValleyMethod,charge);
                isotopeProbeStandard.isotopeNumber = i;
                CgProbe isotopeProbeGreedy = LipidomicsAnalyzer.calculateOneArea(chrom, scan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
                isotopeProbeGreedy.isotopeNumber = i;
                isotopeProbeGreedy.greedyProbe = true;
                double[] isotopeRegionAreas = this.calculateLowerUpperArea(i, theoreticalAreaValue,aProbe.greedyProbe);
                double lowerArea = isotopeRegionAreas[0];
                double upperArea = isotopeRegionAreas[1];
                if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea) || (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                    (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                  CgProbe rightProbe = isotopeProbe;
                  if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                      (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)||
                      (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                    if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)&&
                        isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeGreedy);
                    } else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }
                  }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea){
                    rightProbe = isotopeProbe;
                  }else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                    rightProbe = isotopeProbeStandard;
                  }else if (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                    rightProbe = isotopeProbeGreedy;
                  }
                  this.checkDuplicate(oneisotopeProbes, rightProbe);
                  oneisotopeProbes.add(rightProbe);
                  isOneOfInterest = true;
                }
              }
            }
            if (isOneOfInterest){
              Vector<CgProbe> oneisotopeProbeOK = new Vector<CgProbe>();
              for (CgProbe anProbe: oneisotopeProbes){
                if (anProbe.AreaStatus==CgAreaStatus.OK){
                  oneisotopeProbeOK.add(anProbe);
                }
              }
              if (oneisotopeProbeOK.size()>0){
                isoAreasOfOneProbe.put(i,oneisotopeProbeOK);
                isotopicAreas.put(j, isoAreasOfOneProbe);
              }else{
                isOfInterest[j] = false;
                isotopicAreas.remove(j);
              }
            }else{
              isOfInterest[j] = false;
              isotopicAreas.remove(j);
            }             
          }
        }
      }
    }
    biggestArea = 0;
    highestProbe = null;
    for (int i=0;i!=probes.size();i++){
      if (isOfInterest[i]){
        CgProbe aProbe = probes.get(i);
        if (aProbe.AreaStatus==CgAreaStatus.OK){
          if (aProbe.Area>biggestArea){
            biggestArea = aProbe.Area;
            highestProbe = aProbe;
          }
        }
      }
    }
    for (int i=0;i!=probes.size();i++){
      if (isOfInterest[i]){
        CgProbe aProbe = probes.get(i);
        if (aProbe.AreaStatus==CgAreaStatus.OK){
          if (aProbe.Area>(biggestArea/200)){
            if (this.checkTimeSpace(highestProbe,aProbe)>=30){
              isOfInterest[i] = false;
              isotopicAreas.remove(i);            
            }
          }else{
            isOfInterest[i] = false;
            isotopicAreas.remove(i);
          }
        }
      }
    }


    if (possibleProbabs.size()>probabs.size()){
      for (int i=probabs.size(); i!=possibleProbabs.size();i++){
        double relativeValue = possibleProbabs.get(i)/mainPeak;
        previousIsoChrom = chrom;
        chrom = new LipidomicsChromatogram(readAChromatogram(mz+(i*massToAdd)/(float)charge, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        chrom.GetMaximumAndAverage();
        Vector<CgProbe> allowedPrevIsos = this.getAllowedPrevIsos(isotopicAreas,(i-1),isOfInterest);
        for (int j=0; j!=probes.size();j++){
          if (isOfInterest[j]){
            CgProbe aProbe = probes.get(j);
            Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneProbe = new Hashtable<Integer,Vector<CgProbe>>();
            if (isotopicAreas.containsKey(j))
              isoAreasOfOneProbe = isotopicAreas.get(j);
            int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, chrom);
            result = this.calculateOneAreaWithoutPrevIso(chrom,previousIsoChrom, scan,multiplicationFactor,allowedPrevIsos,false,charge);
            theProbes = (Vector<CgProbe>)result.get(0);
            double theoreticalAreaValue =  aProbe.Area*relativeValue;
            Vector<CgProbe> oneisotopeProbes = new Vector<CgProbe>();
            boolean isOneOfInterest = false;
            for (CgProbe isotopeProbe : theProbes){
              isotopeProbe.isotopeNumber = i;
              if (isotopeProbe.AreaStatus==CgAreaStatus.OK){
                CgProbe isotopeProbeStandard = LipidomicsAnalyzer.calculateOneArea(chrom, scan, LipidomicsDefines.StandardValleyMethod,charge);
                isotopeProbeStandard.isotopeNumber = i;
                CgProbe isotopeProbeGreedy = LipidomicsAnalyzer.calculateOneArea(chrom, scan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
                isotopeProbeGreedy.isotopeNumber = i;
                isotopeProbeGreedy.greedyProbe = true;
                double[] isotopeRegionAreas = this.calculateLowerUpperArea(i, theoreticalAreaValue,aProbe.greedyProbe);
                double lowerArea = isotopeRegionAreas[0];
                double upperArea = isotopeRegionAreas[1];
                if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea) || (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                    (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                  CgProbe rightProbe = isotopeProbe;
                  if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                      (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)||
                      (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                    if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)&&
                        isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeGreedy);
                    } else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }
                  }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea){
                    rightProbe = isotopeProbe;
                  }else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                    rightProbe = isotopeProbeStandard;
                  }else if (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                    rightProbe = isotopeProbeGreedy;
                  }
                  this.checkDuplicate(oneisotopeProbes, rightProbe);
                  oneisotopeProbes.add(rightProbe);
                  isOneOfInterest = true;
                }

              }
            }
            if (isOneOfInterest){
              Vector<CgProbe> oneisotopeProbeOK = new Vector<CgProbe>();
              for (CgProbe anProbe: oneisotopeProbes){
                if (anProbe.AreaStatus==CgAreaStatus.OK){
                  oneisotopeProbeOK.add(anProbe);
                }
              }
              isoAreasOfOneProbe.put(i,oneisotopeProbeOK);
              isotopicAreas.put(j, isoAreasOfOneProbe);
            }

          }
        }  
      }
    }    
    return isotopicAreas;

  }

  @SuppressWarnings("unchecked")
  public Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> processByMzAndProbabs(float mz, int charge,//float mzTolerance, 
      Vector<Double>probabs, Vector<Double>possibleProbabs, int msLevel, boolean negative) throws CgException{
//    return this.processByMzProbabsAndPossibleRetentionTime(mz, //mzTolerance, 
//        -1, -1, -1, -1, probabs, possibleProbabs);
    this.coarseChromMzTolerance_ = LipidomicsConstants.getCoarseChromMzTolerance(mz);
    float massToAdd = LipidomicsConstants.getNeutronMass();
    LipidomicsChromatogram mainChrom = new LipidomicsChromatogram(readAChromatogram(mz, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
    mainChrom.GetMaximumAndAverage();
    LipidomicsChromatogram previousIsoChrom = new LipidomicsChromatogram(readAChromatogram(mz-massToAdd/(float)charge, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
    previousIsoChrom.GetMaximumAndAverage();
    double multiplicationFactor = (probabs.get(0)/probabs.get(1))*(2f/3f);


    int currentScan =0;
    CgProbe probe;
    Vector<CgProbe> allProbes = new Vector<CgProbe>();
    @SuppressWarnings("rawtypes")
    Vector result = null;
    Vector<CgProbe> theProbes = null;
    while (currentScan<mainChrom.Value.length){
      result = this.calculateOneAreaWithoutPrevIso(mainChrom,previousIsoChrom, currentScan,multiplicationFactor,new Vector<CgProbe>(),true,charge);
      theProbes = (Vector<CgProbe>)result.get(0);
      int highestScan = (Integer)result.get(2);
      for (int i=0;i!=theProbes.size();i++){
        probe = theProbes.get(i);
        probe.isotopeNumber = 0;
        if (probe.AreaStatus==CgAreaStatus.OK){
          this.checkDuplicate(allProbes,probe);
          if (!probe.greedyProbe&&probe.AreaStatus==CgAreaStatus.Duplicate){
            probe = LipidomicsAnalyzer.calculateOneArea(mainChrom, currentScan, LipidomicsDefines.StandardValleyMethod,charge);
          }
          int isPossibleIsotope = 0;
          if (probe.greedyProbe)
            isPossibleIsotope = isFromOtherIsotope(probe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>());
          else
            isPossibleIsotope = isFromOtherIsotope(probe, LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>());
          //          if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP){
//            probe = LipidomicsAnalyzer.calculateOneArea(mainChrom, currentScan, LipidomicsDefines.GreedySteepnessReductionMethod);
//            probe.greedyProbe = true;
//          }
          if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP && (!probe.greedyProbe)){
//            System.out.println("1111111");
            CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(mainChrom, currentScan, LipidomicsDefines.GreedySteepnessReductionMethod);
            greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(probe,greedyProbe,mainChrom,previousIsoChrom, multiplicationFactor,new Vector<CgProbe>(),charge);
            //            if (!isGreedyFromOtherIsotope(greedyProbe,probe,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>())){
//            if (isFromOtherIsotope(greedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom, multiplicationFactor,new Vector<CgProbe>())!= LipidomicsAnalyzer.ISOTOPE_PRESENT){
            if (greedyProbe!=null){
            //              System.out.println("22222222222");
              greedyProbe.greedyProbe = true;
              probe = greedyProbe;
            }//else{
//              System.out.println("333333");
//              isPossibleIsotope=LipidomicsAnalyzer.ISOTOPE_PRESENT;
//            }
          }

          //          if (!isFromOtherIsotope(probe, LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,new Vector<CgProbe>())){
          if (isPossibleIsotope!=LipidomicsAnalyzer.ISOTOPE_PRESENT){
            this.checkDuplicate(allProbes,probe);
            allProbes.add(probe);
          }
          Vector<CgProbe> noDuplicateProbes = new Vector<CgProbe>();
          for (CgProbe aProbe: allProbes){
            if (aProbe.AreaStatus == CgAreaStatus.OK)
              noDuplicateProbes.add(aProbe);
          }
          allProbes = noDuplicateProbes;
//          int upperValleyIndex = this.findIndexByTime(probe.UpperValley,mainChrom);
//          if (upperValleyIndex>currentScan)
//            currentScan = upperValleyIndex;
          if (highestScan>currentScan)
            currentScan = highestScan;
        }
      }
      currentScan++;
    }
    // exclude Areas which are probably from a different isotopic pattern
//    LipidomicsChromatogram chrom2 = this.readAChromatogram(mz-NEUTRON_MASS, mzTolerance, mzTolerance);
//    chrom2.GetMaximumAndAverage();
    boolean[] isOfInterest = new boolean[allProbes.size()];
    for (int i=0; i!=allProbes.size();i++){
      isOfInterest[i] = true;
    }
//    for (int j=0; j!=allProbes.size();j++){
//      CgProbe aProbe = allProbes.get(j);
//      int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, chrom2);
//      CgProbe otherIsotopeProbe = LipidomicsAnalyzer.calculateOneArea(chrom2, scan, LipidomicsDefines.EnhancedValleyMethod);
//      CgProbe otherIsotopeProbe2 = LipidomicsAnalyzer.calculateOneArea(chrom2, scan, LipidomicsDefines.StandardValleyMethod);
//      CgProbe otherOfInterest;
//      //This is just a check if too much has been quantified with EnhancedValley
//      if (otherIsotopeProbe.Area>(100*otherIsotopeProbe2.Area)){
//        otherOfInterest=otherIsotopeProbe2;
//      }else{
//        otherOfInterest=otherIsotopeProbe;
//      }
//      if (otherOfInterest.Area>(aProbe.Area*multiplicationFactor)){
//        if (otherOfInterest.Peak>aProbe.LowerValley&&otherOfInterest.Peak<aProbe.UpperValley){
//          isOfInterest[j] = false;
//        }  
//      }
//    }

    
    double mainPeak = 1;
    Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicAreas = new Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>>();
    for (int j=0; j!=allProbes.size();j++){
      if (isOfInterest[j]){
        Vector<CgProbe> areas = new Vector<CgProbe>();
        areas.add(allProbes.get(j));
        Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneIso = new Hashtable<Integer,Vector<CgProbe>>();
        isoAreasOfOneIso.put(0, areas);
        isotopicAreas.put(j, isoAreasOfOneIso);
      }
    }
    for (int i=0;i!=probabs.size();i++){
      previousIsoChrom = mainChrom;
      if (i==0){
        mainPeak = probabs.get(i);
      }else{
        double relativeValue = probabs.get(i)/mainPeak;
        mainChrom = new LipidomicsChromatogram(readAChromatogram(mz+(i*massToAdd)/(float)charge, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        mainChrom.GetMaximumAndAverage();
        Vector<CgProbe> allowedPrevIsos = this.getAllowedPrevIsos(isotopicAreas,(i-1),isOfInterest);
        for (int j=0; j!=allProbes.size();j++){
          if (isOfInterest[j]){
            CgProbe aProbe = allProbes.get(j);
            Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneProbe = new Hashtable<Integer,Vector<CgProbe>>();
            if (isotopicAreas.containsKey(j))
              isoAreasOfOneProbe = isotopicAreas.get(j);
            int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, mainChrom);
            result = this.calculateOneAreaWithoutPrevIso(mainChrom,previousIsoChrom, scan,multiplicationFactor,allowedPrevIsos,false,charge);
            theProbes = (Vector<CgProbe>)result.get(0);
            boolean isOneOfInterest = false;
            double theoreticalAreaValue =  aProbe.Area*relativeValue;
            Vector<CgProbe> oneisotopeProbes = new Vector<CgProbe>();
            for (CgProbe isotopeProbe : theProbes){
              isotopeProbe.isotopeNumber = i;
              if (isotopeProbe.AreaStatus==CgAreaStatus.OK){
                CgProbe isotopeProbeStandard = LipidomicsAnalyzer.calculateOneArea(mainChrom, scan, LipidomicsDefines.StandardValleyMethod,charge);
                isotopeProbeStandard.isotopeNumber = i;
                CgProbe isotopeProbeGreedy = LipidomicsAnalyzer.calculateOneArea(mainChrom, scan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
                isotopeProbeGreedy.isotopeNumber = i;
                isotopeProbeGreedy.greedyProbe = true;
                double[] isotopeRegionAreas = this.calculateLowerUpperArea(i, theoreticalAreaValue,aProbe.greedyProbe);
                double lowerArea = isotopeRegionAreas[0];
                double upperArea = isotopeRegionAreas[1];
                if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea) || (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                    (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                  CgProbe rightProbe = isotopeProbe;
                  if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                      (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)||
                      (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                    if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)&&
                        isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeGreedy);
                    } else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }
                  }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea){
                    rightProbe = isotopeProbe;
                  }else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                    rightProbe = isotopeProbeStandard;
                  }else if (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                    rightProbe = isotopeProbeGreedy;
                  }
                  this.checkDuplicate(oneisotopeProbes, rightProbe);
                  oneisotopeProbes.add(rightProbe);
                  isOneOfInterest = true;
                }
              }
            }
            if (isOneOfInterest){
              Vector<CgProbe> oneisotopeProbeOK = new Vector<CgProbe>();
              for (CgProbe anProbe: oneisotopeProbes){
                if (anProbe.AreaStatus==CgAreaStatus.OK){
                  oneisotopeProbeOK.add(anProbe);
                }
              }
              if (oneisotopeProbeOK.size()>0){
                isoAreasOfOneProbe.put(i,oneisotopeProbeOK);
                isotopicAreas.put(j, isoAreasOfOneProbe);
              }else{
                isOfInterest[j] = false;
                isotopicAreas.remove(j);
              }
            }else{
              isOfInterest[j] = false;
              isotopicAreas.remove(j);
            }             
          }
        }
      }
    }
    float biggestArea = 0;
    CgProbe highestProbe = null;
    for (int i=0;i!=allProbes.size();i++){
      if (isOfInterest[i]){
        CgProbe aProbe = allProbes.get(i);
        if (aProbe.AreaStatus==CgAreaStatus.OK){
          if (aProbe.Area>biggestArea){
            biggestArea = aProbe.Area;
            highestProbe = aProbe;
          }
        }
      }
    }
    for (int i=0;i!=allProbes.size();i++){
      if (isOfInterest[i]){
        CgProbe aProbe = allProbes.get(i);
        if (aProbe.AreaStatus==CgAreaStatus.OK){
          if (aProbe.Area>(biggestArea/200)){
            if (this.checkTimeSpace(highestProbe,aProbe)>=30){
              isOfInterest[i] = false;
              isotopicAreas.remove(i);            
            }
          }else{
            isOfInterest[i] = false;
            isotopicAreas.remove(i);
          }
        }
      }
    }


    if (possibleProbabs.size()>probabs.size()){
      for (int i=probabs.size(); i!=possibleProbabs.size();i++){
        double relativeValue = possibleProbabs.get(i)/mainPeak;
        previousIsoChrom = mainChrom;
        mainChrom = new LipidomicsChromatogram(readAChromatogram(mz+(i*massToAdd)/(float)charge, this.coarseChromMzTolerance_, this.coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        mainChrom.GetMaximumAndAverage();
        Vector<CgProbe> allowedPrevIsos = this.getAllowedPrevIsos(isotopicAreas,(i-1),isOfInterest);
        for (int j=0; j!=allProbes.size();j++){
          if (isOfInterest[j]){
            CgProbe aProbe = allProbes.get(j);
            Hashtable<Integer,Vector<CgProbe>> isoAreasOfOneProbe = new Hashtable<Integer,Vector<CgProbe>>();
            if (isotopicAreas.containsKey(j))
              isoAreasOfOneProbe = isotopicAreas.get(j);
            int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, mainChrom);
            result = this.calculateOneAreaWithoutPrevIso(mainChrom,previousIsoChrom, scan,multiplicationFactor,allowedPrevIsos,false,charge);
            theProbes = (Vector<CgProbe>)result.get(0);
            double theoreticalAreaValue =  aProbe.Area*relativeValue;
            Vector<CgProbe> oneisotopeProbes = new Vector<CgProbe>();
            boolean isOneOfInterest = false;
            for (CgProbe isotopeProbe : theProbes){
              isotopeProbe.isotopeNumber = i;
              if (isotopeProbe.AreaStatus==CgAreaStatus.OK){
                CgProbe isotopeProbeStandard = LipidomicsAnalyzer.calculateOneArea(mainChrom, scan, LipidomicsDefines.StandardValleyMethod,charge);
                isotopeProbeStandard.isotopeNumber = i;
                CgProbe isotopeProbeGreedy = LipidomicsAnalyzer.calculateOneArea(mainChrom, scan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
                isotopeProbeGreedy.isotopeNumber = i;
                isotopeProbeGreedy.greedyProbe = true;
                double[] isotopeRegionAreas = this.calculateLowerUpperArea(i, theoreticalAreaValue,aProbe.greedyProbe);
                double lowerArea = isotopeRegionAreas[0];
                double upperArea = isotopeRegionAreas[1];
                if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea) || (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                    (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                  CgProbe rightProbe = isotopeProbe;
                  if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)||
                      (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)||
                      (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea)){
                    if ((lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea)&&
                        isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeStandard);
                    }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbe,isotopeProbeGreedy);
                    } else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK&&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea && isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&  lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                      rightProbe = this.checkCloserToExpected(aProbe,isotopeProbeStandard,isotopeProbeGreedy);
                    }
                  }else if (lowerArea<isotopeProbe.Area&&isotopeProbe.Area<upperArea){
                    rightProbe = isotopeProbe;
                  }else if (isotopeProbeStandard.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeStandard.Area&&isotopeProbeStandard.Area<upperArea){
                    rightProbe = isotopeProbeStandard;
                  }else if (isotopeProbeGreedy.AreaStatus==CgAreaStatus.OK &&lowerArea<isotopeProbeGreedy.Area&&isotopeProbeGreedy.Area<upperArea){
                    rightProbe = isotopeProbeGreedy;
                  }
                  this.checkDuplicate(oneisotopeProbes, rightProbe);
                  oneisotopeProbes.add(rightProbe);
                  isOneOfInterest = true;
                }

              }
            }
            if (isOneOfInterest){
              Vector<CgProbe> oneisotopeProbeOK = new Vector<CgProbe>();
              for (CgProbe anProbe: oneisotopeProbes){
                if (anProbe.AreaStatus==CgAreaStatus.OK){
                  oneisotopeProbeOK.add(anProbe);
                }
              }
              isoAreasOfOneProbe.put(i,oneisotopeProbeOK);
              isotopicAreas.put(j, isoAreasOfOneProbe);
            }

          }
        }  
      }
    }    
    return isotopicAreas;
  }

  private double[] calculateLowerUpperArea(int isotopeNumber, double theoreticalAreaValue,boolean greedyProbe){
    double[] areas = new double[2];
    double lowerArea = 0d;
    double upperArea = 0d;
/*    if (isotopeNumber==1){
      lowerArea = (2*theoreticalAreaValue)/3;
      upperArea = (4*theoreticalAreaValue)/3;
    }else{*/
      lowerArea = theoreticalAreaValue/3;
      upperArea = 2*theoreticalAreaValue;
//    }
//    if (greedyProbe){
//      lowerArea = theoreticalAreaValue/4;
//      upperArea = 3*theoreticalAreaValue;
//    }
    areas[0] = lowerArea;
    areas[1] = upperArea;
    return areas;
  }

  /** This method checks if there are any peaks from previous chromatograms; if so the peak is neglected;
   *  Furthermore it gives start/stop scans where to continue
   *  The returned probes always start with the one with the lowest retention time*/
  @SuppressWarnings({ "unchecked", "rawtypes" })
  protected Vector calculateOneAreaWithoutPrevIso(LipidomicsChromatogram chrom,LipidomicsChromatogram previousIsoChrom, int mainScan, 
      double multiplicationFactor,Vector<CgProbe> allowedPrevIsos, boolean respectIntCutoff,int charge){
    Vector result = new Vector();
    Vector<CgProbe> probes = new Vector<CgProbe>();
    if (!respectIntCutoff || chrom.Value[mainScan][1]>this.highestIntensity_*this.generalBasePeakCutoff_){
    CgProbe valleyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.EnhancedValleyMethod,charge);
    CgProbe standardProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.StandardValleyMethod,charge);
//    System.out.println(mainScan);

    if (valleyProbe.AreaStatus==CgAreaStatus.OK/*&&standardProbe.AreaStatus==CgAreaStatus.OK*/){
      //check if all three are the same
      if (standardProbe.AreaStatus==CgAreaStatus.OK&&(standardProbe.Area-1)<valleyProbe.Area&&valleyProbe.Area<(standardProbe.Area+1)){
//        System.out.println("1111111");
        int isPossibleIsotope = isFromOtherIsotope(standardProbe,LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,allowedPrevIsos);
//        System.out.println(isPossibleIsotope);
        if (isPossibleIsotope != LipidomicsAnalyzer.NO_ISOTOPE_PRESENT){
          CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
//          System.out.println("greedyProbe.Area: "+greedyProbe.Area);
          greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(standardProbe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,allowedPrevIsos,charge);
          if (greedyProbe!=null){
            valleyProbe = greedyProbe;
            valleyProbe.greedyProbe = true; 
            isPossibleIsotope = LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP;
          }
        }
        if (isPossibleIsotope == LipidomicsAnalyzer.ISOTOPE_PRESENT){
          valleyProbe.AreaStatus = CgAreaStatus.TooSmall;
        }
        probes.add(valleyProbe);
        if (valleyProbe.AreaStatus == CgAreaStatus.OK){
          if (this.highestArea_<valleyProbe.Area)
            this.highestArea_ = valleyProbe.Area;
          if (this.highestIntensity_<valleyProbe.getHighestIntensity())
            this.highestIntensity_ = valleyProbe.getHighestIntensity();         
        }
        result.add(probes);
        result.add(LipidomicsAnalyzer.findIndexByTime(valleyProbe.LowerValley,chrom));
        result.add(LipidomicsAnalyzer.findIndexByTime(valleyProbe.UpperValley,chrom));
      }else{
//        System.out.println("2222222222");
        Vector<CgProbe> allProbes = new Vector<CgProbe>();
        Vector<CgProbe> okStandardProbes = this.calculateProbesWithinSpecificBoundaries(valleyProbe, standardProbe,chrom, LipidomicsDefines.StandardValleyMethod,charge);
        // now check if  one of the probes is from a different isotope
        Vector<Boolean> fromDifferentIsotope = new Vector<Boolean>();
        boolean notOneFromDifferentIsotope = true;
        Vector<CgProbe> okStandardProbesNew = new Vector<CgProbe>();
        for (CgProbe aProbe: okStandardProbes){         
          int isPossibleIsotope = isFromOtherIsotope(aProbe, LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,allowedPrevIsos);
          if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP){
            CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, LipidomicsAnalyzer.findIndexByTime(aProbe.Peak,chrom), LipidomicsDefines.GreedySteepnessReductionMethod,charge);
            greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(standardProbe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,allowedPrevIsos,charge);
            if (greedyProbe!=null/*&&!isGreedyFromOtherIsotope(greedyProbe,aProbe,previousIsoChrom,multiplicationFactor,allowedPrevIsos)*/){
              greedyProbe.greedyProbe = true;
              okStandardProbesNew.add(greedyProbe);
              notOneFromDifferentIsotope = false;
            }else{
              okStandardProbesNew.add(aProbe);
            }
            fromDifferentIsotope.add(false);
          } else if (isPossibleIsotope != LipidomicsAnalyzer.NO_ISOTOPE_PRESENT){
            fromDifferentIsotope.add(true);
            notOneFromDifferentIsotope = false;
            okStandardProbesNew.add(aProbe);
          }else{
            fromDifferentIsotope.add(false);
            okStandardProbesNew.add(aProbe);
          }
          
        }
        okStandardProbes = okStandardProbesNew;
        if (notOneFromDifferentIsotope){
          allProbes.add(valleyProbe);
        }else{
          for (int i=0; i!=okStandardProbes.size();i++){
            if (!fromDifferentIsotope.get(i)){
              int count = this.getSortPosition(allProbes, okStandardProbes.get(i));
              this.checkDuplicate(allProbes,okStandardProbes.get(i));
              allProbes.add(count,okStandardProbes.get(i));
            }else{
              CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, LipidomicsAnalyzer.findIndexByTime(okStandardProbes.get(i).Peak,chrom), LipidomicsDefines.GreedySteepnessReductionMethod,charge);
              greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(standardProbe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,allowedPrevIsos,charge);
              if (greedyProbe!=null/*&&!isGreedyFromOtherIsotope(greedyProbe,okStandardProbes.get(i),previousIsoChrom,multiplicationFactor,allowedPrevIsos)*/){
                int count = this.getSortPosition(allProbes, greedyProbe);
                this.checkDuplicate(allProbes,greedyProbe);
                greedyProbe.greedyProbe = true;
                allProbes.add(count,greedyProbe);
              }
            }
          }
        }
        for (CgProbe aProbe: allProbes){
          if (aProbe.AreaStatus==CgAreaStatus.OK)
            probes.add(aProbe);
        }
        result.add(probes);
        result.add(LipidomicsAnalyzer.findIndexByTime(valleyProbe.LowerValley,chrom));
        result.add(LipidomicsAnalyzer.findIndexByTime(valleyProbe.UpperValley,chrom));    
      }
    }else if (standardProbe.AreaStatus==CgAreaStatus.OK){
//      System.out.println("33333333333333");
      int isPossibleIsotope = isFromOtherIsotope(standardProbe,LipidomicsDefines.StandardValleyMethod,previousIsoChrom,multiplicationFactor,allowedPrevIsos);
      
      if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP){
        CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
        greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(standardProbe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,allowedPrevIsos,charge);
        if (greedyProbe!=null){

//        if (isFromOtherIsotope(greedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom, multiplicationFactor,allowedPrevIsos)!= LipidomicsAnalyzer.ISOTOPE_PRESENT){
//        if (!isGreedyFromOtherIsotope(greedyProbe,standardProbe,previousIsoChrom,multiplicationFactor,allowedPrevIsos)){
          standardProbe = greedyProbe;
          standardProbe.greedyProbe = true;               
        }//else
         // isPossibleIsotope = LipidomicsAnalyzer.ISOTOPE_PRESENT;
      }
      if (isPossibleIsotope == LipidomicsAnalyzer.ISOTOPE_PRESENT){
        standardProbe.AreaStatus = CgAreaStatus.TooSmall;
      }
      probes.add(standardProbe);
      result.add(probes);
      result.add(LipidomicsAnalyzer.findIndexByTime(standardProbe.LowerValley,chrom));
      result.add(LipidomicsAnalyzer.findIndexByTime(standardProbe.UpperValley,chrom));      
    }else{
//      System.out.println("444444444444");
//      if (isFromOtherIsotope(valleyProbe,LipidomicsDefines.EnhancedValleyMethod,previousIsoChrom,multiplicationFactor,allowedPrevIsos)){
//        valleyProbe.AreaStatus = CgAreaStatus.TooSmall;
//      }
      int isPossibleIsotope = isFromOtherIsotope(standardProbe,LipidomicsDefines.EnhancedValleyMethod,previousIsoChrom,multiplicationFactor,allowedPrevIsos);
      if (isPossibleIsotope == LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP){
        CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
        greedyProbe = this.selectCorrespondingGreedyProbeIsGreedy(standardProbe,greedyProbe,chrom,previousIsoChrom, multiplicationFactor,allowedPrevIsos,charge);
        if (greedyProbe!=null){

//        if (!isGreedyFromOtherIsotope(greedyProbe,valleyProbe,previousIsoChrom,multiplicationFactor,allowedPrevIsos)){
//        if (isFromOtherIsotope(greedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom, multiplicationFactor,allowedPrevIsos)!= LipidomicsAnalyzer.ISOTOPE_PRESENT){
          valleyProbe = greedyProbe;
          valleyProbe.greedyProbe = true;               
        }//else
         // isPossibleIsotope = LipidomicsAnalyzer.ISOTOPE_PRESENT;
      }
      if (isPossibleIsotope == LipidomicsAnalyzer.ISOTOPE_PRESENT){
        valleyProbe.AreaStatus = CgAreaStatus.TooSmall;
      }
      probes.add(valleyProbe);
      result.add(probes);
      result.add(LipidomicsAnalyzer.findIndexByTime(valleyProbe.LowerValley,chrom));
//      System.out.println(LipidomicsAnalyzer.findIndexByTime(valleyProbe.UpperValley,chrom));
      result.add(LipidomicsAnalyzer.findIndexByTime(valleyProbe.UpperValley,chrom));  
    }
    }else{
      result.add(probes);
      result.add(mainScan);
//      System.out.println(LipidomicsAnalyzer.findIndexByTime(valleyProbe.UpperValley,chrom));
      result.add(mainScan); 
    }
    return result;
  }
  
  private CgProbe selectCorrespondingGreedyProbeIsGreedy(CgProbe standardProbe, CgProbe greedyProbe, LipidomicsChromatogram chrom, LipidomicsChromatogram previousIsoChrom, double multiplicationFactor,Vector<CgProbe> allowedPrevIsos,int charge){
    Vector<CgProbe> okStandardProbes = this.calculateProbesWithinSpecificBoundaries(standardProbe, greedyProbe,chrom, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
    Vector<CgProbe> probesWithoutIso = new Vector<CgProbe>();
    for (CgProbe newGreedyProbe : okStandardProbes){
      if (newGreedyProbe.Area>(standardProbe.Area/10)){
        if (!isGreedyFromOtherIsotope(newGreedyProbe, standardProbe, previousIsoChrom,multiplicationFactor,allowedPrevIsos)){
          probesWithoutIso.add(newGreedyProbe);
        }
      }
    }
    CgProbe returnProbe = null;
    for (CgProbe probeWithoutIso : probesWithoutIso){
      if (returnProbe==null){
        returnProbe = probeWithoutIso;
      }else{
        if (probeWithoutIso.Area>returnProbe.Area)
          returnProbe = probeWithoutIso;
      }       
    }
    return returnProbe;
  }
  
//  private CgProbe selectCorrespondingGreedyProbe(CgProbe standardProbe, CgProbe greedyProbe, LipidomicsChromatogram chrom, LipidomicsChromatogram previousIsoChrom, double multiplicationFactor,Vector<CgProbe> allowedPrevIsos){
//    Vector<CgProbe> okStandardProbes = this.calculateProbesWithinSpecificBoundaries(standardProbe, greedyProbe,chrom, LipidomicsDefines.GreedySteepnessReductionMethod);
//    Vector<CgProbe> probesWithoutIso = new Vector<CgProbe>();
//    for (CgProbe newGreedyProbe : okStandardProbes){
//      if (isFromOtherIsotope(newGreedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod,previousIsoChrom, multiplicationFactor,allowedPrevIsos)!= LipidomicsAnalyzer.ISOTOPE_PRESENT){
//        probesWithoutIso.add(newGreedyProbe);
//      }
//    }
//    CgProbe returnProbe = null;
//    for (CgProbe probeWithoutIso : probesWithoutIso){
//      if (returnProbe==null){
//        returnProbe = probeWithoutIso;
//      }else{
//        if (probeWithoutIso.Area>returnProbe.Area)
//          returnProbe = probeWithoutIso;
//      }       
//    }
//    return returnProbe;
//  }
  
  private Vector<CgProbe> calculateProbesWithinSpecificBoundaries(CgProbe valleyProbe, CgProbe standardProbe,LipidomicsChromatogram chrom, int mode, int charge){
    int lowerIndex = LipidomicsAnalyzer.findIndexByTime(valleyProbe.LowerValley,chrom);
    int upperIndex = LipidomicsAnalyzer.findIndexByTime(valleyProbe.UpperValley,chrom);
    int currentScan = LipidomicsAnalyzer.findIndexByTime(standardProbe.LowerValley,chrom);
    Vector<CgProbe> allStandardProbes = new Vector<CgProbe>();
    if (standardProbe.AreaStatus==CgAreaStatus.OK){
      allStandardProbes.add(standardProbe);
    }
    CgProbe currentValidProbe = standardProbe;
    CgProbe probe;
    while (currentScan>=lowerIndex){
      probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, mode,charge);
      if (mode==LipidomicsDefines.GreedySteepnessReductionMethod)
        probe.greedyProbe = true;
      if (probe.AreaStatus==CgAreaStatus.OK){
        this.checkDuplicate(allStandardProbes,probe);
        if (currentValidProbe==null||probe.AreaStatus==CgAreaStatus.OK){
          currentValidProbe = probe;
        }
        allStandardProbes.add(probe);
        int lowerValleyIndex = LipidomicsAnalyzer.findIndexByTime(probe.LowerValley,chrom);
        if (lowerValleyIndex<currentScan)
          currentScan = lowerValleyIndex;
      }
      currentScan--;
    }
    currentScan = LipidomicsAnalyzer.findIndexByTime(standardProbe.UpperValley,chrom);
    while (currentScan<=upperIndex){
      probe = LipidomicsAnalyzer.calculateOneArea(chrom, currentScan, mode,charge);
      if (mode==LipidomicsDefines.GreedySteepnessReductionMethod)
        probe.greedyProbe = true;
      if (probe.AreaStatus==CgAreaStatus.OK){
        this.checkDuplicate(allStandardProbes,probe);
        if (currentValidProbe==null||probe.AreaStatus==CgAreaStatus.OK){
          currentValidProbe = probe;
        }
        allStandardProbes.add(probe);
        int upperValleyIndex = LipidomicsAnalyzer.findIndexByTime(probe.UpperValley,chrom);
        if (upperValleyIndex>currentScan)
          currentScan = upperValleyIndex;
      }
      currentScan++;
    }
    Vector<CgProbe> okStandardProbes = new Vector<CgProbe>();
    for (CgProbe aProbe : allStandardProbes){
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        okStandardProbes.add(aProbe);
      }  
    }
    return okStandardProbes;
  }

  private int getSortPosition(Vector<CgProbe> allProbes,CgProbe okStandardProbe){
    //for sorting
    int count = 0;
    for (int j=0;j!=allProbes.size();j++){
      if (okStandardProbe.LowerValley<allProbes.get(j).LowerValley){
        break;
      }else{
        count = j;
      }
    }
    return count;
  }
  
  private boolean isGreedyFromOtherIsotope(CgProbe greedyProbe, CgProbe standardProbe,LipidomicsChromatogram previousIsoChrom, double multiplicationFactor,Vector<CgProbe> allowedPrevIsos){
    boolean isFromOtherIso = true;
    float prevDistanceGreedy = greedyProbe.Peak-greedyProbe.LowerValley;
    float postDistanceGreedy = greedyProbe.UpperValley-greedyProbe.Peak;
    float prevDistanceStandard = standardProbe.Peak-standardProbe.LowerValley;
    float postDistanceStandard = standardProbe.UpperValley-standardProbe.Peak;
    boolean checkPrevIsoGreedy = false;
    // this checks if there is 10% times more distance between the peak and the valley in greedy mode
    // if there is a difference then an overlap is possible
    if (prevDistanceStandard>(1.1*prevDistanceGreedy))
      checkPrevIsoGreedy = true;
    if (postDistanceStandard>(1.1*postDistanceGreedy))
      checkPrevIsoGreedy = true;
 //   System.out.println(checkPrevIsoGreedy);
    if (checkPrevIsoGreedy){
      int isPossibleIsotope = this.isFromOtherIsotope(greedyProbe, LipidomicsDefines.GreedySteepnessReductionMethod, previousIsoChrom, multiplicationFactor, allowedPrevIsos);
//      System.out.println("isPossibleIsotope: "+isPossibleIsotope);
      if (isPossibleIsotope != LipidomicsAnalyzer.ISOTOPE_PRESENT)
        isFromOtherIso = false;
    }
    
    return isFromOtherIso;
  }
  
//  private Vector<CgProbe> checkPeakInBetween(Vector<CgProbe> okProbes,LipidomicsChromatogram previousIsoChrom,double multiplicationFactor,Vector<CgProbe> allowedPrevIsos){
//    System.out.println("I am checking the probes against one another");
//    Vector<CgProbe> checkedProbes = new Vector<CgProbe>();
//    Hashtable<Integer,Boolean> isOk = new Hashtable<Integer,Boolean>();
//    for (int i=0; i!=okProbes.size();i++){
//      System.out.println(i);
//      if (!isOk.containsKey(i)){
//        isOk.put(i, true);
//      }
//      for (int j=(i+1); j!=okProbes.size();j++){
//        CgProbe probe1 = okProbes.get(i);
//        CgProbe probe2 = okProbes.get(j);
//        System.out.println(i+";"+j);
//        if (probe1.greedyProbe||probe2.greedyProbe){          
//          CgProbe otherIsotopeProbe1 = Analyzer.calculateOneArea(previousIsoChrom, Analyzer.findIndexByTime(probe1.Peak, previousIsoChrom), LipidomicsDefines.StandardValleyMethod);
//          CgProbe otherIsotopeProbe2 = Analyzer.calculateOneArea(previousIsoChrom, Analyzer.findIndexByTime(probe2.Peak, previousIsoChrom), LipidomicsDefines.StandardValleyMethod);
//          System.out.println(otherIsotopeProbe1.Peak+";"+otherIsotopeProbe1.Area);
//          System.out.println(otherIsotopeProbe2.Peak+";"+otherIsotopeProbe2.Area);
//          //the isotopic peak is the same and it is between the two peaks
//          if ((otherIsotopeProbe1.Peak-0.1)<otherIsotopeProbe2.Peak&&otherIsotopeProbe2.Peak<(otherIsotopeProbe1.Peak+0.1)&&
//            ((probe1.Peak<otherIsotopeProbe1.Peak&&probe2.Peak>otherIsotopeProbe1.Peak) ||(probe1.Peak>otherIsotopeProbe1.Peak&&probe2.Peak<otherIsotopeProbe1.Peak))){
//            System.out.println("In between");
//            // now I know the isotopic peak is between them
//            // then I have to check if that is a possible isotope value
//            if (otherIsotopeProbe1.Area>(probe1.Area*multiplicationFactor)||otherIsotopeProbe1.Area>(probe2.Area*multiplicationFactor)){
//              if (this.isAllowedPrevIso(otherIsotopeProbe1, allowedPrevIsos)){
//                System.out.println("It is a prev iso and in between");
//                float distance1 = otherIsotopeProbe1.Peak-probe1.Peak;
//                if (distance1<0)
//                  distance1 = distance1*(-1);
//                float distance2 = otherIsotopeProbe1.Peak-probe2.Peak;
//                if (distance2<0)
//                  distance2 = distance2*(-1);
//                if (distance1==distance2){
//                  if (probe1.Peak>probe2.Peak)
//                    isOk.put(j,false);
//                  else
//                    isOk.put(i,false);
//                }else if (distance1>distance2)
//                  isOk.put(j,false);
//                else
//                  isOk.put(i,false);
//              }
//            }
//          }
//        }
//      }
//    }
//    for (int i=0; i!=okProbes.size();i++){
//      if (isOk.get(i))
//        checkedProbes.add(okProbes.get(i));
//    }
//    return checkedProbes;
//  }
  
  private boolean isAllowedPrevIso(CgProbe otherIsotopeProbeInterest, Vector<CgProbe> allowedPrevIsos){
    boolean isAnAllowedIsotope = false;
    for (CgProbe allowedIso : allowedPrevIsos){
      //
//      System.out.println(allowedIso.LowerValley+";"+allowedIso.UpperValley);
//      System.out.println(otherIsotopeProbeInterest.LowerValley+";"+otherIsotopeProbeInterest.UpperValley+";"+allowedIso.greedyProbe);
      if (((allowedIso.LowerValley-1)<otherIsotopeProbeInterest.LowerValley&&otherIsotopeProbeInterest.UpperValley<(allowedIso.UpperValley+1))||
          (allowedIso.greedyProbe&&allowedIso.isCoveredByThisProbe(otherIsotopeProbeInterest))){
        isAnAllowedIsotope = true;
      }
    }
    return isAnAllowedIsotope;
  }
  
  private int isFromOtherIsotope(CgProbe aProbe, int mode,LipidomicsChromatogram previousIsoChrom, double multiplicationFactor,Vector<CgProbe> allowedPrevIsos){
   int hasIsotope = LipidomicsAnalyzer.NO_ISOTOPE_PRESENT;
   int scan = LipidomicsAnalyzer.findIndexByTime(aProbe.Peak, previousIsoChrom);
   CgProbe otherIsotopeProbeInterest;
   if (mode==LipidomicsDefines.StandardValleyMethod || (mode==LipidomicsDefines.GreedySteepnessReductionMethod)){
     otherIsotopeProbeInterest = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, scan, LipidomicsDefines.StandardValleyMethod);
   }else{
     CgProbe otherIsotopeProbe = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, scan, LipidomicsDefines.EnhancedValleyMethod);
     CgProbe otherIsotopeProbe2 = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, scan, LipidomicsDefines.StandardValleyMethod);
   
     //This is just a check if too much has been quantified with EnhancedValley
     if (otherIsotopeProbe.Area>(100*otherIsotopeProbe2.Area)){
       otherIsotopeProbeInterest=otherIsotopeProbe2;
     }else{
       otherIsotopeProbeInterest=otherIsotopeProbe;
     }
   }
   int possibleIsotopeValue = LipidomicsAnalyzer.NO_ISOTOPE_PRESENT;
   if (otherIsotopeProbeInterest.Area>(aProbe.Area*multiplicationFactor)){
     possibleIsotopeValue = LipidomicsAnalyzer.ISOTOPE_PRESENT;
//     if (mode==LipidomicsDefines.GreedySteepnessReductionMethod){
//       System.out.println("Greedy check");
//       System.out.println(otherIsotopeProbeInterest.LowerValley+";"+otherIsotopeProbeInterest.Peak+";"+otherIsotopeProbeInterest.UpperValley);
//       System.out.println(aProbe.LowerValley+";"+aProbe.Peak+";"+aProbe.UpperValley);
//     }
   } else if (otherIsotopeProbeInterest.Area>(aProbe.Area/3)){
//     System.out.println("Possible overlap!!!!!!");
     possibleIsotopeValue = LipidomicsAnalyzer.POSSIBLE_ISOTOPE_OVERLAP;
   }
   if (possibleIsotopeValue!=LipidomicsAnalyzer.NO_ISOTOPE_PRESENT){
//     if (mode==LipidomicsDefines.GreedySteepnessReductionMethod){
//       System.out.println(possibleIsotopeValue);
//       System.out.println(aProbe.Area+";"+otherIsotopeProbeInterest.Peak+";"+aProbe.LowerValley+";"+aProbe.UpperValley+";"+otherIsotopeProbeInterest.LowerValley+";"+otherIsotopeProbeInterest.UpperValley+";"+aProbe.isCoveredByThisProbe(otherIsotopeProbeInterest));
//     }
     if (/*(mode==LipidomicsDefines.GreedySteepnessReductionMethod&&(*/(otherIsotopeProbeInterest.Peak<aProbe.LowerValley||otherIsotopeProbeInterest.Peak>aProbe.UpperValley)&&
         (mode!=LipidomicsDefines.GreedySteepnessReductionMethod||(!aProbe.isCoveredByThisProbe(otherIsotopeProbeInterest)))/*))||
         (otherIsotopeProbeInterest.LowerValley<aProbe.LowerValley&&otherIsotopeProbeInterest.UpperValley<aProbe.LowerValley)||
         (otherIsotopeProbeInterest.UpperValley>aProbe.UpperValley&&otherIsotopeProbeInterest.LowerValley>aProbe.UpperValley)*/);
//     boolean otherIsoPossible = true;
//     if ((otherIsotopeProbeInterest.Peak<aProbe.LowerValley||otherIsotopeProbeInterest.Peak>aProbe.UpperValley)/*)||
//         (mode==LipidomicsDefines.GreedySteepnessReductionMethod&&(!aProbe.isCoveredByThisProbe(otherIsotopeProbeInterest))))||
//         (otherIsotopeProbeInterest.LowerValley<aProbe.LowerValley&&otherIsotopeProbeInterest.UpperValley<aProbe.LowerValley)||
//         (otherIsotopeProbeInterest.UpperValley>aProbe.UpperValley&&otherIsotopeProbeInterest.LowerValley>aProbe.UpperValley)*/){
//       otherIsoPossible = false;
//     }else{
//       if (!this.checkInterPeakDistanceSmall(otherIsotopeProbeInterest,aProbe)&&(mode!=LipidomicsDefines.GreedySteepnessReductionMethod||(!aProbe.isCoveredByThisProbe(otherIsotopeProbeInterest))))
//         otherIsoPossible = false;
//     }
//     if (otherIsoPossible){
     else if (mode==LipidomicsDefines.GreedySteepnessReductionMethod&&(!aProbe.greedyFragment)&&this.isInOuterThird(otherIsotopeProbeInterest, aProbe)){
//       if (aProbe.Area>100000000){
//         System.out.println("Hello I am checking a greedy probe: "+aProbe.Area);
//         System.out.println(aProbe.LowerValley+";"+aProbe.UpperValley);
//         System.out.println(aProbe.greedyFragment+";"+this.checkInterPeakDistanceSmall(otherIsotopeProbeInterest,aProbe)+(!aProbe.isCoveredByThisProbe(otherIsotopeProbeInterest)));
//         &&(!aProbe.greedyFragment)&&(!this.checkInterPeakDistanceSmall(otherIsotopeProbeInterest,aProbe)));
//       }
     }
     else{
       boolean isAnAllowedIsotope = this.isAllowedPrevIso(otherIsotopeProbeInterest, allowedPrevIsos);
       if (!isAnAllowedIsotope)
         hasIsotope = possibleIsotopeValue;
     }  
   }
   return hasIsotope;
  }
  
  /** checks the how close the borders of the enhanced valley peak and the standard valley peak are to the one
   * of the 0 isotope; if there is one of the 1.1 times higher than the other, the closer one is taken, if not the
   * preference is given to the enhanced peak
   * @param originalProbe 0 isotope probe
   * @param isotopeProbe other isotopic peak calculated by enhance valley method
   * @param isotopeProbeStandard other isotopic peak calculated by standard valley method
   * @return
   */
  private CgProbe checkCloserToExpected(CgProbe originalProbe, CgProbe isotopeProbe, CgProbe isotopeProbeStandard){
    float totalDistanceEnh =  this.positiveDifference(originalProbe.LowerValley, isotopeProbe.LowerValley)+this.positiveDifference(originalProbe.UpperValley, isotopeProbe.UpperValley);    
    float totalDistanceStand =  this.positiveDifference(originalProbe.LowerValley, isotopeProbeStandard.LowerValley)+this.positiveDifference(originalProbe.UpperValley, isotopeProbeStandard.UpperValley);    
    float ratio = totalDistanceEnh/totalDistanceStand;
    if (ratio<(1/1.1)||1.1<ratio){
      if (totalDistanceEnh<totalDistanceStand){
        return isotopeProbe;
      }else{
        return isotopeProbeStandard;
      }
    }else{
      return isotopeProbe;
    }
  }
  
  private CgProbe checkCloserToExpected(CgProbe originalProbe, CgProbe isotopeProbe, CgProbe isotopeProbeStandard, CgProbe greedyProbe){
    CgProbe cgProbe1 = this.checkCloserToExpected(originalProbe, isotopeProbe, isotopeProbeStandard);
    CgProbe cgProbe2 = this.checkCloserToExpected(originalProbe, isotopeProbeStandard, greedyProbe);
    return this.checkCloserToExpected(originalProbe, cgProbe1, cgProbe2);
  }

  
  private float positiveDifference(float nr1, float nr2){
    float distance = nr1-nr2;
    if (distance>0){
      return distance;
    }else{
      return (-1)*distance;
    }
  }
  
  private Vector<CgProbe> getAllowedPrevIsos(Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isoAreasOfOneProbe,int isotope,boolean[] isOfInterest){
    Vector<CgProbe> allowedPrevIsos = new Vector<CgProbe>();
    for (int i=0;i!=isOfInterest.length;i++){
      if (isOfInterest[i]){
        Hashtable<Integer,Vector<CgProbe>> areaIsoHashOfOneProbe = isoAreasOfOneProbe.get(i);
        if (areaIsoHashOfOneProbe.containsKey(isotope)){
          allowedPrevIsos.addAll(areaIsoHashOfOneProbe.get(isotope));
        }  
      }  
    }  
    return allowedPrevIsos;
  }
  
  private float checkTimeSpace(CgProbe highestProbe,CgProbe aProbe){
    float timeSpace=0;
    if (!(aProbe.Peak-0.1<highestProbe.Peak&&highestProbe.Peak<aProbe.Peak+0.1)){
      if (aProbe.Peak>highestProbe.Peak){
        timeSpace=aProbe.LowerValley-highestProbe.UpperValley;
      }else{
        timeSpace=highestProbe.LowerValley-aProbe.UpperValley;
      }
    }
    return timeSpace;    
  }
  
  private boolean isInOuterThird(CgProbe otherIsotopeProbeInterest, CgProbe aProbe){
    boolean isInOutherThird = false;
    float peakDistance = otherIsotopeProbeInterest.Peak-aProbe.Peak;
    if (peakDistance<0)
      peakDistance = peakDistance*-1;
    if (otherIsotopeProbeInterest.Peak<aProbe.Peak){
      float outerRegionThreshold = aProbe.LowerValley+(aProbe.Peak-aProbe.LowerValley)/3;
      if (otherIsotopeProbeInterest.Peak<=outerRegionThreshold)
        isInOutherThird = true;
    }else{
      float outerRegionThreshold = aProbe.UpperValley-(aProbe.UpperValley-aProbe.Peak)/3;
      if (otherIsotopeProbeInterest.Peak>=outerRegionThreshold)
        isInOutherThird = true;
    }
        return isInOutherThird;
  }

  /** this method finds a 3D-peak for a specific scan*/
  public Probe3D detectPeakThreeD(LipidomicsChromatogram coarseChrom, int mainScan, boolean respectIntCutoff, int charge, int msLevel) throws CgException,QuantificationException{
    Probe3D probe3D = null;
    if (!respectIntCutoff || (coarseChrom.Value[mainScan][1]>(this.highestIntensity_*this.generalBasePeakCutoff_) && coarseChrom.Value[mainScan][2]>this.lowerIntensityThreshold_)){
    // First a peak (2D) from the coarse chromatogram is calculated
    CgProbe coarseProbe = LipidomicsAnalyzer.calculateOneArea(coarseChrom, mainScan, LipidomicsDefines.GreedySteepnessImproved,charge);
    checkIfSmoothingDeviatedPeakPosition(coarseChrom,coarseProbe);
//    System.out.println("Coarse Probe: "+coarseProbe.LowerValley+" ; "+coarseProbe.Peak+" ; "+coarseProbe.UpperValley+";"+coarseChrom.getLoSteepnessPoint()+";"+coarseChrom.getUpSteepnessPoint()+";"+mainScan);
    // if the coarse probe is OK continue
    if (coarseProbe.AreaStatus==CgAreaStatus.OK/*&&(!useSameCgHashFor3D_||!inSameCgHash(coarseProbe))*/){
//      System.out.println("coarseProbe.Area: "+coarseProbe.Area);
      // Second an m/z profile at the timepoint of the coarse-probe peak is extracted 
      CgProbe rightProfileProbe = this.getProfileProbe(coarseProbe,coarseChrom.Mz,charge,msLevel);
//      System.out.println(rightProfileProbe.LowerValley+";"+rightProfileProbe.Peak+";"+rightProfileProbe.UpperValley+";"+rightProfileProbe.AreaStatus);
      if (rightProfileProbe.AreaStatus==CgAreaStatus.OK){
        // in the profile probe the correct m/z value for the peak is detected; with this 
        // value a 2D-chromatogram with with a broad and a small mz Range is extracted
        @SuppressWarnings("rawtypes")
        Vector theProbes = this.getSmallAndBroadProbe(rightProfileProbe, mainScan, coarseChrom, charge,msLevel);
        CgProbe finalProbe = (CgProbe)theProbes.get(1);
//        System.out.println("Final Probe: "+finalProbe.LowerValley+" ; "+finalProbe.Peak+" ; "+finalProbe.UpperValley);
        // the peak has not been found where the original one was; hence a small overlapping one
        // is assumed; recalculate with the new one and shift the peak
        // the condition to assume this is that the peak of the previous coarse probe is not within a specific range
        // after some repeats the procedure is stopped
        int count = 0;
        while ((coarseProbe.Peak<(finalProbe.LowerValley+finalProbeTimeCompTolerance_)||(finalProbe.UpperValley-finalProbeTimeCompTolerance_)<coarseProbe.Peak ||
            finalProbe.Mz<(rightProfileProbe.LowerValley+finalProbeMzCompTolerance_)||finalProbe.Mz>(rightProfileProbe.UpperValley-finalProbeMzCompTolerance_))&& count<5){
          coarseProbe = finalProbe;
          rightProfileProbe = this.getProfileProbe(finalProbe,coarseChrom.Mz,charge,msLevel);
          theProbes = this.getSmallAndBroadProbe(rightProfileProbe, mainScan, (LipidomicsChromatogram)theProbes.get(2), charge,msLevel);
          finalProbe =(CgProbe)theProbes.get(1);
          count++;
        }
//        System.out.println(coarseProbe.Peak+";"+finalProbe.LowerValley+";"+finalProbe.UpperValley);
//        System.out.println(finalProbe.Mz+";"+rightProfileProbe.LowerValley+";"+rightProfileProbe.UpperValley);
        // if the conditions of above are still not fulfilled after some repeats the peak is assumed to be not valid 
        if (count==5&&(coarseProbe.Peak<(finalProbe.LowerValley+finalProbeTimeCompTolerance_)||(finalProbe.UpperValley-finalProbeTimeCompTolerance_)<coarseProbe.Peak ||
            finalProbe.Mz<(rightProfileProbe.LowerValley+finalProbeMzCompTolerance_)||finalProbe.Mz>(rightProfileProbe.UpperValley-finalProbeMzCompTolerance_))){
          coarseProbe.AreaStatus=CgAreaStatus.TooSmall;
          probe3D = new Probe3D(coarseProbe);
          //System.out.println("4444444444444");
          return probe3D;
        }
        //System.out.println("rightProfileProbe : "+finalProbe.LowerValley+";"+finalProbe.UpperValley);
        //System.out.println("finalProbe : "+rightProfileProbe.LowerValley+";"+rightProfileProbe.UpperValley+" ; "+finalProbe.Area+" ; ");
        //System.out.println(rightProfileProbe.AreaStatus+";"+finalProbe.AreaStatus);
        
        // this is just to introduce a fixed maximum length for the m/z direction
        // the reason is that a peak can sometimes creep on the ground without a big change of the steepness
        // or intensity; this prevents the creeping profile peaks
        if (rightProfileProbe.LowerValley<(finalProbe.Mz-profileMaxRange_)){
          rightProfileProbe.LowerValley = finalProbe.Mz-profileMaxRange_;
        }  
        if (rightProfileProbe.UpperValley>(finalProbe.Mz+profileMaxRange_))
          rightProfileProbe.UpperValley = finalProbe.Mz+profileMaxRange_;
        if (rightProfileProbe.AreaStatus==CgAreaStatus.OK&&finalProbe.AreaStatus==CgAreaStatus.OK){
          // now that everything seems to be OK we can calculate an ellipse which form the ground of the 3D peak
          float[] ellipse = Calculator.calculateEllipseFrom4Points(coarseProbe.Peak, rightProfileProbe.LowerValley,
            coarseProbe.Peak, rightProfileProbe.UpperValley, 
          finalProbe.LowerValley, finalProbe.Mz, finalProbe.UpperValley, finalProbe.Mz);
          // if the highest intensity of the peak is not within the peak borders, the probe is set to the coarse probe peak
          if (finalProbe.Peak<(finalProbe.LowerValley+finalProbeTimeCompTolerance_)||(finalProbe.UpperValley-finalProbeTimeCompTolerance_)<finalProbe.Peak)
            finalProbe.Peak = coarseProbe.Peak;
//          System.out.println("x0: "+ellipse[0]+"; y: "+ellipse[1]+"; a: "+ellipse[2]+"; b: "+ellipse[3]);
//          if (Float.isNaN(ellipse[3])||(!Float.isInfinite(ellipse[3]))){
//            System.out.println("!!!!!!!!!!!!!!!!!!!!!!!");
//            ellipse[3] = (rightProfileProbe.UpperValley-rightProfileProbe.LowerValley)/2;
//          }
          if (ellipse[3]>profileMaxRange_)
            ellipse[3] = profileMaxRange_;
          // now a new object containing the 3D-parameters and the information about a possible overlap is generated
          probe3D =this.getOverlapPosition(finalProbe, rightProfileProbe, coarseProbe, ellipse);
          // now a chromatogram containing just the intensities within the ellips borders is extracted
          LipidomicsChromatogram forQuant = new LipidomicsChromatogram(readJustIntensitiesOfInterest(probe3D,msLevel));
          
          forQuant.Background = finalProbe.Background;    
          // just the raw intensities are extracted and stored as area; no smoothing is required
          forQuant.getAreaRaw();
          probe3D.Area = forQuant.Area;
          if (probe3D.Area>highestArea_)
            highestArea_ = probe3D.Area;
          if (forQuant.getHighestIntensity()>highestIntensity_)
            highestIntensity_ = forQuant.getHighestIntensity();
          
//        if (probe3D.Area>10000000f)
//        System.out.println(probe3D.Area+";"+probe3D.isOverlapBefore()+";"+probe3D.isOverlapLater()+";"+coarseChrom.Value[mainScan][0]);
//        System.out.println(coarseProbe.LowerValley+";"+coarseProbe.Peak+";"+coarseProbe.UpperValley);
//        System.out.println("ReturnedProbe: "+probe3D.LowerValley+";"+probe3D.Peak+";"+probe3D.UpperValley);
//        System.out.println(probe3D.Mz+";"+probe3D.LowerMzBand+";"+probe3D.UpperMzBand);
//        if (useSameCgHashFor3D_&& coarseProbe.LowerValley<probe3D.Peak&&probe3D.Peak<coarseProbe.UpperValley)
////        if (useSameCgHashFor3D_ && (probe3D.LowerValley+(probe3D.Peak-probe3D.LowerValley)/3)<coarseProbe.Peak&&
////            coarseProbe.Peak<(probe3D.UpperValley-(probe3D.UpperValley-probe3D.Peak)/3)){
////          this.sameCgHash_.add(coarseProbe);
////          this.same3DProbes_.add(probe3D);
////        }  
//        System.out.println("probe3D.Area: "+probe3D.Area);
//        System.out.println("probe3D.Area: "+probe3D.Area);
/*      XYSeries series1 = new XYSeries("Chrom small");
      XYSeries series2 = new XYSeries("Chrom broad 5");
      XYSeries series3 = new XYSeries("Chrom broad 10");
      for (int j=0;j!=cx1.ScanCount;j++){
        series1.add(cx1.Value[j][0],cx1.Value[j][2]);
//        series2.add(cx2.Value[j][0],cx2.Value[j][2]);
        series3.add(cx3.Value[j][0],cx3.Value[j][2]);
      }
      XYSeriesCollection coll1 = new XYSeriesCollection(series1);
      XYSeriesCollection coll2 = new XYSeriesCollection(series2);
      XYSeriesCollection coll3 = new XYSeriesCollection(series3);
      JFreeChart chart1 = ChartFactory.createXYLineChart("Chrom small","time","Intensity",coll1,PlotOrientation.VERTICAL,true,false,false);
      JFreeChart chart2 = ChartFactory.createXYLineChart("Chrom broad 5","time","Intensity",coll2,PlotOrientation.VERTICAL,true,false,false);
      JFreeChart chart3 = ChartFactory.createXYLineChart("Chrom broad 10","time","Intensity",coll3,PlotOrientation.VERTICAL,true,false,false);
      BufferedOutputStream stream;
      try {
        stream = new BufferedOutputStream(new FileOutputStream("F:\\lipidomics\\testProfiles\\"+3+"_chrom_1.JPG"));
        ImageIO.write(chart1.createBufferedImage(1000, 700), "JPEG", stream);
        stream.close();
        stream = new BufferedOutputStream(new FileOutputStream("F:\\lipidomics\\testProfiles\\"+3+"_chrom_2.JPG"));
        ImageIO.write(chart2.createBufferedImage(1000, 700), "JPEG", stream);
        stream.close();
        stream = new BufferedOutputStream(new FileOutputStream("F:\\lipidomics\\testProfiles\\"+3+"_chrom_3.JPG"));
        ImageIO.write(chart3.createBufferedImage(1000, 700), "JPEG", stream);
        stream.close();
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      XYSeries series1 = new XYSeries("Chrom small");
      for (int j=0;j!=forQuant.ScanCount;j++){
        series1.add(forQuant.Value[j][0],forQuant.Value[j][1]);
      }
      XYSeriesCollection coll1 = new XYSeriesCollection(series1);
      JFreeChart chart1 = ChartFactory.createXYLineChart("Chrom 3D extract","time","Intensity",coll1,PlotOrientation.VERTICAL,true,false,false);
      BufferedOutputStream stream;
      try {
        stream = new BufferedOutputStream(new FileOutputStream("F:\\lipidomics\\testProfiles\\3D_chrom_1.JPG"));
        ImageIO.write(chart1.createBufferedImage(1000, 700), "JPEG", stream);
        stream.close();
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }*/
//    }
        // if the last profile probe or the final probe is not OK, cast the coarse probe to a 3D probe and say the area is too small   
        }else{
//          System.out.println("33333333333333333");
          if (rightProfileProbe.AreaStatus!=CgAreaStatus.OK)
            coarseProbe.AreaStatus = rightProfileProbe.AreaStatus;
          if (finalProbe.AreaStatus!=CgAreaStatus.OK)
            coarseProbe.AreaStatus = finalProbe.AreaStatus;
          probe3D = new Probe3D(coarseProbe);
        }
     // if the profile probe is not OK, cast the coarse probe to a 3D probe and say the area is too small   
      }else{
//        System.out.println("222222222222222222");
        coarseProbe.AreaStatus = rightProfileProbe.AreaStatus;
        probe3D = new Probe3D(coarseProbe);
      }
    // if the coarse probe is not OK, cast it to a 3D probe and say the area is too small  
    }else{
//      System.out.println("111111111111111111");
//      if (coarseProbe.AreaStatus!=CgAreaStatus.OK){
        coarseProbe.AreaStatus = CgAreaStatus.TooSmall;
        probe3D = new Probe3D(coarseProbe);
//      }else if (inSameCgHash(coarseProbe)){
//        probe3D = new Probe3D(this.getCorresponding3DProbe(coarseProbe));
//        probe3D.AreaStatus = CgAreaStatus.AllreadyUsed;
//      }
    }
    }else{
      probe3D = new Probe3D(new CgProbe(0,charge));
      probe3D.AreaStatus = CgAreaStatus.TooSmall;
    }
    return probe3D;
  }
  
  protected Probe3D getOverlapPosition(CgProbe finalProbe, CgProbe profileProbe, CgProbe oldChromProbe, float[] ellipse){
    boolean before = false;
    boolean later = false;
    boolean lessMz = false;
    boolean higherMz = false;
    // overlaps in the m/z direction are currently not taken into account since it would cause a massive time consumption
    if (profileProbe.hasLoOverlap())
      lessMz = true;
    if (profileProbe.hasUpOverlap())
      higherMz = true;
    // an overlap is assumed of one of the probes hase definitely an overlap or if the distance between
    // the peak and the border deviates heavily between the old and the final probe
    if (oldChromProbe.hasLoOverlap()||finalProbe.hasLoOverlap())
      before = true;
    if (((oldChromProbe.Peak-oldChromProbe.LowerValley)/(finalProbe.Peak-finalProbe.LowerValley))>overlapDistanceDeviationFactor_)
      before = true;
    if (oldChromProbe.hasUpOverlap()||finalProbe.hasUpOverlap())
      later = true;
    if (((oldChromProbe.UpperValley-oldChromProbe.Peak)/(finalProbe.UpperValley-finalProbe.Peak)>overlapDistanceDeviationFactor_))
      later = true;
    Probe3D probe3D = new Probe3D(finalProbe,profileProbe,lessMz,higherMz,before,later,ellipse);
    return probe3D;
  }
  

  /** returns an m/z profile for a specific chromatogram peak*/
  private CgProbe getProfileProbe(CgProbe aProbe, float originalMz, int charge, int msLevel) throws CgException{
    // before the profile is extracted from the file it is checked if such a profile
    // probe has already been calculated; if yes, it can be assumed that it will lead to 
    // the same profile (since the time where the peak has been found is the same);
    // hence the profiles are stored in a hash; if it is same coarse probe, just the profile
    // of the hash is taken and returned
    int probePositionInHash = this.getProbePositionInHash(aProbe);
    if (this.useSameCgHashFor3D_&&probePositionInHash!=-1){
//      System.out.println("Getting from Hash");
      return this.fromProbeToProfile_.get(probePositionInHash);
    }else{
//      System.out.println("Reading from file");
      LipidomicsChromatogram  profile = readASingleProfile(aProbe,msLevel,true);
      
//      this.printChromaToFile(profile, "D:\\Christer\\20170310\\test2\\profsmall.JPG");
//      this.printChromaToFile(profile, "D:\\Christer\\20170310\\test2\\profraw.JPG",1);
      
      profile.GetMaximumAndAverage();
      // these are specific paramaters for the Greedy-Algorithm for detection of peak borders
      profile.greedySteepnessChange1_ = profileSteepnessChange1_;
      profile.greedyIntensityCutoff1_ = profileIntensityCutoff1_;
      profile.greedySteepnessChange2_ = profileSteepnessChange2_;
      profile.greedyIntensityCutoff2_ = profileIntensityCutoff2_;
      // for the peak detection cut-off value relative to the highest intensity can be defined
      profile.intensityCutoff_ = profileGeneralIntCutoff_;
      // it could happen that the peak you are looking for is a little bit shifted in m/z direction
      // therefore this method is called which takes the first peak that is OK in the m/z neighborhood
      int scanStart = profile.Value.length/2;
      float dist = Math.abs(profile.Value[scanStart][0]-originalMz);
      if ((scanStart-1)>-1 && Math.abs(profile.Value[scanStart-1][0]-originalMz)<dist){
        while ((scanStart-1)>-1 && Math.abs(profile.Value[scanStart-1][0]-originalMz)<dist){
          scanStart--;
          dist = Math.abs(profile.Value[scanStart][0]-originalMz);
        }
      } else if ((scanStart+1)<profile.Value.length && Math.abs(profile.Value[scanStart+1][0]-originalMz)<dist){
        while ((scanStart+1)<profile.Value.length && Math.abs(profile.Value[scanStart+1][0]-originalMz)<dist){
          scanStart++;
          dist = Math.abs(profile.Value[scanStart][0]-originalMz);
        }

      }
      profile.findPeakInNeighbourhood(scanStart, 1, LipidomicsDefines.GreedySteepnessImproved,50);
    
////    profile.smoothingBorderCorrection(3);
////    profile.FindPeak(profile.Value.length/2, 1, LipidomicsDefines.StandardValleyMethod);
////    profile.FindPeak(profile.Value.length/2, 1, LipidomicsDefines.EnhancedValleyMethod);
      CgProbe rightProfileProbe = new CgProbe(0,charge);
      rightProfileProbe = LipidomicsAnalyzer.copyResults(rightProfileProbe, profile, true);
//    System.out.println(rightProfileProbe.Peak+";"+(originalMz-(mzRange*2)/3)+";"+originalMz+";"+(originalMz+((mzRange*2)/3)));
      // if the peak is to far out of the
      rightProfileProbe.Peak = rightProfileProbe.Peak+profileSmoothingCorrection_;
      //System.out.println(originalMz+";"+rightProfileProbe.Peak+";"+profilePeakAcceptanceRange_);
      //System.out.println("DeltaMz: "+(rightProfileProbe.Peak-originalMz));
      
      if ((rightProfileProbe.Peak)<(originalMz-profilePeakAcceptanceRange_)||(originalMz+profilePeakAcceptanceRange_)<(rightProfileProbe.Peak))
        rightProfileProbe.AreaStatus = CgAreaStatus.TooSmall;
      else if ((rightProfileProbe.UpperValley-rightProfileProbe.LowerValley)<profileMzMinRange_)
        rightProfileProbe.AreaStatus = CgAreaStatus.MzRangeTooShort;
      if (this.useSameCgHashFor3D_){
        this.fromProbeToProfile_.put(this.chromProbes_.size(), rightProfileProbe);
        this.chromProbes_.add(aProbe);
      }
      //System.out.println("Borders "+rightProfileProbe.LowerValley+";"+rightProfileProbe.Peak+";"+rightProfileProbe.UpperValley+";"+rightProfileProbe.hasUpOverlap());
      return rightProfileProbe;
    }
  }

  /** reads an m/z profile from the chrom file and smooths it*/
  protected Vector<CgChromatogram> readProfiles(Vector<CgProbe> probes, float mzTolerance, float timeTolerance, float maxTimeDeviation,
      float mzSmoothRange, int smoothRepeats, float meanSmoothRange, int meanSmoothRepeats, int msLevel, boolean smooth) throws CgException{
    LipidomicsChromReader lReader = (LipidomicsChromReader) reader_;
    Vector<CgChromatogram> chroms =  lReader.readProfiles(probes, mzTolerance, timeTolerance, maxTimeDeviation,0,0, msLevel, sav_gol_jni_);
    if (!smooth) return chroms;
    for (CgChromatogram chrom : chroms){
      LipidomicsChromatogram lChrom = new LipidomicsChromatogram(chrom);
      boolean copyRawData = true;
      if (meanSmoothRange>0){
        lChrom.smoothMean(meanSmoothRange, meanSmoothRepeats,copyRawData);
        copyRawData = false;
      }
      if( useCuda_ ){
        lChrom.Smooth(mzSmoothRange, smoothRepeats, copyRawData, sav_gol_jni_);
      } else {
        lChrom.Smooth(mzSmoothRange, smoothRepeats, copyRawData);
      }
    }
    return chroms;
  }
   
  @SuppressWarnings({ "unchecked", "rawtypes" })
  private Vector getSmallAndBroadProbe(CgProbe rightProfileProbe, int mainScan, LipidomicsChromatogram coarseChrom, int charge, int msLevel)throws CgException,QuantificationException{
    //System.out.println("rightProfileProbe.Peak: "+rightProfileProbe.Peak);
    LipidomicsChromatogram smallChroma = null;
    LipidomicsChromatogram broadChroma = null;
    // if the same profile probe has been detected, the result will be the 
    // extraction of the same small and broad chromatogram; hence a hashtable
    // is used to save time
    int profilePositionInHash = this.getProfilePositionInHash(rightProfileProbe);
    if (this.useSameCgHashFor3D_&&profilePositionInHash!=-1){
//      System.out.println("Chrom from Hash");
      smallChroma = this.fromProfileToSmallChrom_.get(profilePositionInHash);
      broadChroma = this.fromProfileToBroadChrom_.get(profilePositionInHash);
    }else{
//      System.out.println("Chrom from File");
      // this section reads a small and a broad chromatogram for the analysis
      // the m/z-range for the small chromatogram is settable, while the range
      // of the broad one corresponds to the borders calculated from the profile
      float smallChromSmoothRange = smallChromSmoothRange_*60f;
      float broadChromSmoothRange = broadChromSmoothRange_*60f;
      smallChroma = new LipidomicsChromatogram(readAChromatogram(rightProfileProbe.Peak, smallChromMzRange_, smallChromMzRange_,msLevel,chromSmoothRange_,smallChromSmoothRepeats_,chromSmoothRange_,smallChromMeanSmoothRepeats_,rightProfileProbe.Mz-smallChromSmoothRange,rightProfileProbe.Mz+smallChromSmoothRange));
      smallChroma.GetMaximumAndAverage();
      smallChroma.intensityCutoff_ = smallChromIntensityCutoff_;
      //broadChroma = new LipidomicsChromatogram(readAChromatogram(rightProfileProbe.Peak, rightProfileProbe.Peak-rightProfileProbe.LowerValley, rightProfileProbe.UpperValley-rightProfileProbe.Peak,msLevel,chromSmoothRange_,broadChromSmoothRepeats_,rightProfileProbe.Mz-broadChromSmoothRange,rightProfileProbe.Mz+broadChromSmoothRange));
      broadChroma = new LipidomicsChromatogram(readAChromatogram(rightProfileProbe.Peak, rightProfileProbe.Peak-rightProfileProbe.LowerValley, rightProfileProbe.UpperValley-rightProfileProbe.Peak,msLevel,chromSmoothRange_,broadChromSmoothRepeats_,chromSmoothRange_,broadChromMeanSmoothRepeats_,rightProfileProbe.Mz-broadChromSmoothRange,rightProfileProbe.Mz+broadChromSmoothRange));
      broadChroma.GetMaximumAndAverage();
      broadChroma.intensityCutoff_ = broadChromIntensityCutoff_;
      // this puts the chromatograms on the hash for further usage
      if (this.useSameCgHashFor3D_){
        this.fromProfileToSmallChrom_.put(this.profileProbes_.size(), smallChroma);
        this.fromProfileToBroadChrom_.put(this.profileProbes_.size(), broadChroma);
        this.profileProbes_.add(rightProfileProbe);
      }
    }
////    smallChroma.FindPeak(mainScan, 1, LipidomicsDefines.GreedySteepnessImproved);
    smallChroma.FindPeak(mainScan, 1, LipidomicsDefines.GreedySteepnessImproved);
    smallChroma.GetBackground(50);
    smallChroma.GetAreaAndTime();
//  this.printChromaToFile(profile, "E:\\lipidomics\\20130502\\12_1_profsmall.JPG");
//  this.printChromaToFile(profile, "E:\\lipidomics\\20130502\\12_1_profraw.JPG",1);

//    this.printChromaToFile(smallChroma, "D:\\Experiment1\\QTOF\\smallChromraw.JPG",1);
//    this.printChromaToFile(smallChroma, "D:\\Experiment1\\QTOF\\smallChrom.JPG");
//    
//    this.printChromaToFile(broadChroma, "D:\\Experiment1\\QTOF\\broadChromraw.JPG",1);
//    this.printChromaToFile(broadChroma, "D:\\Experiment1\\QTOF\\broadChrom.JPG");
    if (!smallChroma.Good) throw new QuantificationException("The exact chromatogram (small) does not show satisfying intensities");
    CgProbe smallProbe = new CgProbe(0,charge);
    smallProbe = LipidomicsAnalyzer.copyResults(smallProbe, smallChroma, false);
//    System.out.println("smallProbe: "+smallProbe.LowerValley/60f+" ; "+smallProbe.Peak/60f+" ; "+smallProbe.UpperValley/60f+" ; "+smallProbe.AreaStatus);

    
//    if (smallProbe.AreaStatus==CgAreaStatus.OK){
    broadChroma.FindPeak(mainScan, 1, LipidomicsDefines.GreedySteepnessPointImproved);
//    System.out.println(broadChroma.Value[broadChroma.LoValley][0]/60f+" ; "+broadChroma.Value[broadChroma.UpValley][0]/60f);
    broadChroma.Peak = smallChroma.Peak;
    
    // in the next to paragraphes is checked if the lower valley or the upper valley of the small valley
    // should be taken as lower and upper valley or not
    // the criteria are if there is a steepness point is found or if the distance difference is big
    // then an overlap must be assumed
    boolean takeLoVallyOfSmall = true;
    if (/*coarseChrom.getLoSteepnessPoint()!=-1&&*/(coarseChrom.LoValley>smallChroma.LoValley)){
      float distanceSmall = smallChroma.Peak-smallChroma.LoValley;
      float distanceBig = coarseChrom.Peak-coarseChrom.LoValley;
//      System.out.println(distanceSmall+";"+distanceBig);
      if ((coarseChrom.getLoSteepnessPoint()!=-1&&distanceSmall>1.1*distanceBig)||distanceSmall>1.5*distanceBig){
        takeLoVallyOfSmall = false;
//        System.out.println("Possible overlap in lower region");
      }
    }
//    if (takeLoVallyOfSmall)
//    if (broadChroma.LoValley<smallChroma.LoValley)
      broadChroma.LoValley = smallChroma.LoValley;
    
      
//    else
//      broadChroma.LoValley = coarseChrom.LoValley;
  
    boolean takeUpVallyOfSmall = true;
    if (/*coarseChrom.getUpSteepnessPoint()!=-1&&*/(coarseChrom.UpValley<smallChroma.UpValley)){
      float distanceSmall = smallChroma.UpValley-smallChroma.Peak;
      float distanceBig = coarseChrom.UpValley-coarseChrom.Peak;
      if ((coarseChrom.getUpSteepnessPoint()!=-1&&distanceSmall>1.1*distanceBig)||distanceSmall>1.5*distanceBig){
        takeUpVallyOfSmall = false;
//        System.out.println("Possible overlap in upper region");
      }         
    }
//    if (broadChroma.UpValley>smallChroma.UpValley)
      broadChroma.UpValley = smallChroma.UpValley;
    
    
    broadChroma.GetBackground(50);
    broadChroma.GetAreaAndTime();
    if (!broadChroma.Good) throw new QuantificationException("The exact chromatogram (broad)  does not show satisfying intensities");
    CgProbe finalProbe = new CgProbe(0,charge);

    finalProbe = LipidomicsAnalyzer.copyResults(finalProbe, broadChroma, false);
//    System.out.println("finalProbe: "+finalProbe.LowerValley+" ; "+finalProbe.Peak+" ; "+finalProbe.UpperValley);
//    System.out.println(takeLoVallyOfSmall+";"+takeUpVallyOfSmall);
//    System.out.println(takeLoVallyOfSmall+";"+takeUpVallyOfSmall+";"+smallProbe.LowerValley+";"+smallProbe.Peak+";"+smallProbe.UpperValley);
    
    // if on one side is an overlap and the peak borders of the small chromatogram should not be taken
    if ((!takeLoVallyOfSmall)||(!takeUpVallyOfSmall)){
      broadChroma.greedySteepnessChange1_ = broadChromSteepnessChangeNoSmall_;
//      System.out.println("I am finding now the peak for the broad chrom");
//      broadChroma.FindPeak(mainScan, 1, LipidomicsDefines.GreedySteepnessImproved);
////      broadChroma.FindPeak(smallChroma.Peak, 1, LipidomicsDefines.GreedySteepnessImproved);
//      System.out.println("I am finding now the peak for the broad chrom finished");
      
//      System.out.println(broadChroma.Value[broadChroma.LoValley][0]+";"+broadChroma.Value[broadChroma.Peak][0]+";"+broadChroma.Value[broadChroma.UpValley][0]);
      int refinedPosition = -1;
      // the lower valley should not be taken
      if (!takeLoVallyOfSmall){
        // find the peak and the borders
        broadChroma.intensityCutoff_ = broadChromIntensityCutoffNoSmall_;
        broadChroma.FindPeak(smallChroma.Peak, 1, LipidomicsDefines.GreedySteepnessImproved);
        // if the lower valley of the broadChroma is higher (peak more narrow) an overlap must be assumed
        // the lower valley of the final probe is now set to the one of the newly calculated one of the broadChroma
        // the boolean loOverlap is set to true
        if (broadChroma.Value[broadChroma.LoValley][0]>finalProbe.LowerValley){
          finalProbe.LowerValley = broadChroma.Value[broadChroma.LoValley][0];
          finalProbe.setLoOverlap(true);
          // the peak position of the small chromatogram could now be out of the newly set
          // range border, if so take the one of the broad chroma because
          // the highest peak intensity must not be outside the peak borders
          if (finalProbe.Peak<finalProbe.LowerValley){
            finalProbe.Peak = broadChroma.Value[broadChroma.Peak][0];
            refinedPosition = broadChroma.Peak;
          }
        }
        // here the whole story is repeated for the main scan  becaus the smallChroma.Peak is sometimes not so reliable
        broadChroma.FindPeak(mainScan, 1, LipidomicsDefines.GreedySteepnessImproved);
        if (broadChroma.Value[broadChroma.LoValley][0]>finalProbe.LowerValley){
          finalProbe.LowerValley = broadChroma.Value[broadChroma.LoValley][0];
          finalProbe.setLoOverlap(true);
          if (finalProbe.Peak<finalProbe.LowerValley){
            finalProbe.Peak = broadChroma.Value[broadChroma.Peak][0];
            refinedPosition = broadChroma.Peak;
          }  
        }
      }
//      System.out.println(broadChroma.Value[broadChroma.UpValley][0]+";"+smallChroma.Value[broadChroma.UpValley][0]);
   // the upper valley should not be taken
      if (!takeUpVallyOfSmall){
     // find the peak and the borders
        broadChroma.intensityCutoff_ = broadChromIntensityCutoffNoSmall_;
        broadChroma.FindPeak(smallChroma.Peak, 1, LipidomicsDefines.GreedySteepnessImproved);
        // if the upper valley of the broadChroma is smaller (peak more narrow) an overlap must be assumed
        // the upper valley of the final probe is now set to the one of the newly calculated one of the broadChroma
        // the boolean upOverlap is set to true
        if (broadChroma.Value[broadChroma.UpValley][0]<finalProbe.UpperValley){
          finalProbe.UpperValley = broadChroma.Value[broadChroma.UpValley][0];
//        finalProbe.UpperValley = (broadChroma.Value[broadChroma.UpValley][0]+smallChroma.Value[smallChroma.UpValley][0])/2;
          finalProbe.setUpOverlap(true);
          // the peak position of the small chromatogram could now be out of the newly set
          // range border, if so take the one of the broad chroma because
          // the highest peak intensity must not be outside the peak borders
          if (finalProbe.Peak>finalProbe.UpperValley){
            finalProbe.Peak = broadChroma.Value[broadChroma.Peak][0];
            refinedPosition = broadChroma.Peak;
          }  
        }
        broadChroma.FindPeak(mainScan, 1, LipidomicsDefines.GreedySteepnessImproved);
     // here the whole story is repeated for the main scan  becaus the smallChroma.Peak is sometimes not so reliable
        if (broadChroma.Value[broadChroma.UpValley][0]<finalProbe.UpperValley){
          finalProbe.UpperValley = broadChroma.Value[broadChroma.UpValley][0];
//        finalProbe.UpperValley = (broadChroma.Value[broadChroma.UpValley][0]+smallChroma.Value[smallChroma.UpValley][0])/2;
          finalProbe.setUpOverlap(true);
          if (finalProbe.Peak>finalProbe.UpperValley){
            finalProbe.Peak = broadChroma.Value[broadChroma.Peak][0];
            refinedPosition = broadChroma.Peak;
          }  
        }
      }
      if (refinedPosition == -1)
        broadChroma.Peak = smallChroma.Peak;
      else
        broadChroma.Peak = refinedPosition;
    }
    Vector theProbes = new Vector();
//    System.out.println(takeLoVallyOfSmall+";"+takeUpVallyOfSmall);
//    System.out.println(finalProbe.LowerValley+";"+finalProbe.UpperValley);
//    this.printChromaToFile(smallChroma, "F:\\lipidomics\\testProfiles\\52_5_0_small_chrom.JPG");
//    this.printChromaToFile(broadChroma, "F:\\lipidomics\\testProfiles\\52_5_0_broad_chrom.JPG");

    // the calculated 2D peaks of the chromatograms are added to the result vector and the latest broadChroma as well
    theProbes.add(smallProbe);
    theProbes.add(finalProbe);
    theProbes.add(broadChroma);
    return theProbes;
  }
  /** this method calculates a the parameters of a peak and checks if the probe is overlapping
   *  to a probe from a previous chromatogram to specify if the peak originates from a different
   *  isotopic distribution*/
  @SuppressWarnings("unused")
  protected Probe3D calculateAreaWithoutPrevIso(LipidomicsChromatogram chrom,LipidomicsChromatogram previousIsoChrom, int mainScan, 
      double multiplicationFactor, boolean respectIntCutoff,int charge,int msLevel) throws CgException{
    // calculates a 3D-probe for the scan
    try{
    Probe3D probe = this.detectPeakThreeD(chrom,mainScan,respectIntCutoff,charge,msLevel);
    if (probe.AreaStatus==CgAreaStatus.OK/*&&(!useSameCgHashFor3D_||!is3DProbeAlreadyUsed(probe))*/){
//      same3DProbes_.add(probe);
      //System.out.println(probe.Area+";"+probe.isOverlapBefore()+";"+probe.isOverlapLater()+";"+probe.isOverlapLowerMz()+";"+probe.isOverlapHigherMz());
      boolean lookLeftSide = true;
      //if there is an overlap on the left side, there is no use to search for potential isotopic peaks on the left
      if (probe.isOverlapBefore())
        lookLeftSide = false;
    //if there is an overlap on the right side, there is no use to search for potential isotopic peaks on the right
      boolean lookRightSide = true;
      if (probe.isOverlapLater())
        lookRightSide = false;
//// considerMz is not used since the 3D peak calculation takes too much time
////      boolean considerMz = false;
      // if there is an overlap in one m/z direction then the isotope m/z value is just allowed on this side
      // if there is no or overlap on both sides this information has no use
      // if there is no overlap in the time axis I do not consider m/z
    ////      if (((probe.isOverlapLowerMz()&&!probe.isOverlapHigherMz())||(!probe.isOverlapLowerMz()&&probe.isOverlapHigherMz()))&&
    ////          (!lookLeftSide||!lookRightSide))
    ////        considerMz = true;
      // collect the possible near Isotopic probes
      Vector<CgProbe> isoProbes = new Vector<CgProbe>();
      CgProbe isoProbe = null;
    ////      if (considerMz){
    ////        isoProbe = this.detectPeakThreeD(previousIsoChrom, mainScan);
    ////      }else{
        isoProbe = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, mainScan, LipidomicsDefines.StandardValleyMethod,charge);
      ////    }
//      System.out.println("MainIsoProbe: "+isoProbe.Area+";"+isoProbe.AreaStatus);
      CgProbe currentValidProbe = null;
      int currentScan = mainScan;
      if (isoProbe.AreaStatus==CgAreaStatus.OK){
        isoProbes.add(isoProbe);
      }
//      if (lookLeftSide){
      currentScan = mainScan-1;
      float valueToFindIndex = probe.LowerValley;
      // if we should look for a peak on the left side we go some steps further to the left
      if (lookLeftSide)
        valueToFindIndex = probe.LowerValley-5;
      int lowerScan = LipidomicsAnalyzer.findIndexByTime(valueToFindIndex,previousIsoChrom);
      if (currentValidProbe!=null){
        int lowerValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.LowerValley,previousIsoChrom);
        if (lowerValleyIndex<currentScan)
          currentScan = lowerValleyIndex;
        currentScan--;
      }
      // no look for isotopic probes on the left side
      while (currentScan>lowerScan){
      ////        if (considerMz)
      ////          isoProbe = this.detectPeakThreeD(previousIsoChrom, mainScan);
      ////       else
          isoProbe = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, currentScan, LipidomicsDefines.StandardValleyMethod,charge);
//          System.out.println("LowerIsoProbe: "+isoProbe.Area+";"+isoProbe.AreaStatus);
        if (isoProbe.AreaStatus==CgAreaStatus.OK){
          isoProbes.add(isoProbe);
          currentValidProbe = isoProbe;
        }
        if (currentValidProbe!=null){
          int lowerValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.LowerValley,previousIsoChrom);
          if (lowerValleyIndex<currentScan)
            currentScan = lowerValleyIndex;
        }
        currentScan-=scanStep_;
      }
      
      currentScan = mainScan+1;
      valueToFindIndex = probe.UpperValley;
   // if we should look for a peak on the right side we go some steps further to the left
      if (lookRightSide)
        valueToFindIndex = probe.UpperValley+5;
      int upperScan = LipidomicsAnalyzer.findIndexByTime(valueToFindIndex,previousIsoChrom);
      if (currentValidProbe!=null){
        int upperValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.UpperValley,previousIsoChrom);
        if (upperValleyIndex>currentScan)
          currentScan = upperValleyIndex;
        currentScan++;
      }
   // no look for isotopic probes on the right side
      while (currentScan<upperScan){
      ////        if (considerMz)
      ////        isoProbe = this.detectPeakThreeD(previousIsoChrom, mainScan);
      ////      else
          isoProbe = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, currentScan, LipidomicsDefines.StandardValleyMethod,charge);
//          System.out.println("UpperIsoProbe: "+isoProbe.Area+";"+isoProbe.AreaStatus);
        if (isoProbe.AreaStatus==CgAreaStatus.OK){
          isoProbes.add(isoProbe);
          currentValidProbe = isoProbe;
        }
        if (currentValidProbe!=null){
          int upperValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.UpperValley,previousIsoChrom);
          if (upperValleyIndex>currentScan)
            currentScan = upperValleyIndex;
        }
        currentScan+=scanStep_;
      }
      // now we have accumulated isotopic probes; we have to check now if they could really cause such a peak
      boolean isFromOtherIso = false;
      for (CgProbe anIsoProbe:isoProbes){
        // the isotopic peak has to have a sufficient intensity (otherwise peaks could be excluded by noise which forms a peak)
        if (anIsoProbe.Area>(probe.Area*multiplicationFactor)){
          // for an overlap one of the peaks top has to be in the range of the other peak
          if ((anIsoProbe.LowerValley<probe.Peak&&probe.Peak<anIsoProbe.UpperValley)||
              probe.LowerValley<anIsoProbe.Peak&&anIsoProbe.Peak<probe.UpperValley){
//            System.out.println(probe.Area+"; 1111111111111111"+"; "+previousIsoChrom.Value[LipidomicsAnalyzer.findIndexByTime(probe.Peak,previousIsoChrom)][2]/previousIsoChrom.Value[LipidomicsAnalyzer.findIndexByTime(anIsoProbe.Peak,previousIsoChrom)][2]);
            // checks if the intensity at the position of the top of the current peak in the previous isotopic chromatogram is in high intensity region
            // something like this has to be done since peak center is not so reliable
            // if the intensity is very low it can not be regarded as overlap
            if (previousIsoChrom.Value[LipidomicsAnalyzer.findIndexByTime(probe.Peak,previousIsoChrom)][2]>(previousIsoChrom.Value[LipidomicsAnalyzer.findIndexByTime(anIsoProbe.Peak,previousIsoChrom)][2]*overlapPossibleIntensityThreshold_)){
              // these are used to check if there is an isotop probe in between
              this.prevIsoProbes_.put(probe.Area, anIsoProbe);
              if (previousIsoChrom.Value[LipidomicsAnalyzer.findIndexByTime(probe.Peak,previousIsoChrom)][2]>(previousIsoChrom.Value[LipidomicsAnalyzer.findIndexByTime(anIsoProbe.Peak,previousIsoChrom)][2]*overlapSureIntensityThreshold_)){
                isFromOtherIso = true;
                break;
              }
              // The probe.isCoveredByThisProbe was too strict!
//              if (this.isOneProbeInOtherInnerThird(probe, anIsoProbe)){ //probe.isCoveredByThisProbe(anIsoProbe)){
//                isFromOtherIso = true;
//                break;
//              }
//              System.out.println(probe.Area+"; 222222222222222222");
              // Since the intensity gave no answer if it is another isotopic peak, 
              // we check now if the peak is nearby the other peak
              // two threshold values are calculated since the  highest peak intensities 
              // could lie sometimes very decentralised, this concerns especially small peaks
              // that overla with big ones
              if (anIsoProbe.Peak<probe.Peak){
//                System.out.println(probe.Area+"; 3333333333333333");
                if (probe.isOverlapBefore()){
//                  System.out.println("xxxxxxxxxxxxx");
                  float lowerThreshold1 = probe.Peak-(probe.Peak-probe.LowerValley)/overlapPeakDistanceDivisor_;
                  float lowerThreshold2 = probe.Peak-(probe.UpperValley-probe.LowerValley)/overlapFullDistanceDivisor_;
                  if (anIsoProbe.Peak>lowerThreshold1||anIsoProbe.Peak>lowerThreshold2){                    
                    isFromOtherIso = true;
                    break;
                  }  
                }else{
                  isFromOtherIso = true;
                  break;
                }  
              }else{
//                System.out.println(probe.Area+"; 444444444444444");
                if (probe.isOverlapLater()){
//                  System.out.println(probe.Area+"; 555555555555555");
                  float upperThreshold1 = probe.Peak+(probe.UpperValley-probe.Peak)/overlapPeakDistanceDivisor_;
                  float upperThreshold2 = probe.Peak+(probe.UpperValley-probe.LowerValley)/overlapFullDistanceDivisor_;
//                  System.out.println(probe.Peak);
//                  float lowerThreshold = anIsoProbe.Peak-(anIsoProbe.Peak-anIsoProbe.LowerValley)/3;
//                  System.out.println(anIsoProbe.Peak+";"+upperThreshold1);
//                  System.out.println(anIsoProbe.Peak+";"+upperThreshold2);
                  if (/*lowerThreshold<probe.Peak||*/anIsoProbe.Peak<upperThreshold1||anIsoProbe.Peak<upperThreshold2){
                    isFromOtherIso = true;
                    break;
                  }  
                }else{
                  isFromOtherIso = true;
                  break;
                }  
              }
            }
          }
        }
      }
//      if (probe.Area>1000000)
//        System.out.println(probe.Area+";"+isFromOtherIso);
      if (isFromOtherIso){
        // if the probe comes from a different isotope this is indicated by a flag
        probe.AreaStatus = CgAreaStatus.OtherIso;
      }
//      }
      // I am not sure if this could lead to problems
    }else{
//      if (probe.Area>0&&probe.Peak>0){
//        CgProbe isoProbe = LipidomicsAnalyzer.calculateOneArea(previousIsoChrom, mainScan, LipidomicsDefines.StandardValleyMethod);
//        if (isoProbe.Area>(probe.Area*multiplicationFactor)&&this.isProbe1InProbe2Third(probe, isoProbe))
//          isoProbe.AreaStatus = CgAreaStatus.OtherIso;
//      }
    }
      return probe;
    
    }catch(QuantificationException qex){
      return null;
    }
  }
  
  
//  private void printChromaToFile(LipidomicsChromatogram cx,String pathName){     
//    this.printChromaToFile(cx, pathName, 2);
//  }
//
//  private void printChromaToFile(LipidomicsChromatogram cx,String pathName,int value){ 
//  
//    XYSeries series1 = new XYSeries("Chrom small");
//    for (int j=0;j!=cx.ScanCount;j++){
////      if (cx.Value[j][0]<782.55f)
//      series1.add(cx.Value[j][0],cx.Value[j][value]);
//    }
//    XYSeriesCollection coll1 = new XYSeriesCollection(series1);
//    JFreeChart chart1 = ChartFactory.createXYLineChart("Chrom small","time","Intensity",coll1,PlotOrientation.VERTICAL,true,false,false);
//    BufferedOutputStream stream;
//    try {
//      stream = new BufferedOutputStream(new FileOutputStream(pathName));
//      ImageIO.write(chart1.createBufferedImage(1000, 700), "JPEG", stream);
//      stream.close();
//    }
//    catch (IOException e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
//    }
//  }

  
  
  /** checks if the peak has already been calculated and if so
   *  returns an identifier to get the chromatogram (-1 corresponds to "not found")
   */ 
  private int getProbePositionInHash(CgProbe chromProbe){
    for (int i=0; i!= this.chromProbes_.size();i++){
      CgProbe probe =  this.chromProbes_.get(i);
      if (this.sameCgProbe(probe, chromProbe))
        return i;
    }
    return -1;
  }
  
  private int getProfilePositionInHash(CgProbe profileProbe){
    for (int i=0; i!= this.profileProbes_.size();i++){
      CgProbe probe =  this.profileProbes_.get(i);
      if (this.sameProfileProbe(probe, profileProbe)){
        return i;
      }  
    }
    return -1;
  }
  
//  private Probe3D getCorresponding3DProbe(CgProbe coarseProbe){
//    for (int i=0; i!= this.sameCgHash_.size();i++){
//      CgProbe probe =  this.sameCgHash_.get(i);
//      if (this.sameCgProbe(probe, coarseProbe))
//        return this.same3DProbes_.get(i);
//    }
//    return null;
//  }
//  
//  private boolean inSameCgHash(CgProbe coarseProbe){
//    boolean inSameHash = false;
//    int count = 0;
//    for (CgProbe probe: this.sameCgHash_){
//      System.out.println("Comparing"+count+": "+coarseProbe.Area+";"+probe.Area+";"+
//          coarseProbe.AreaError+";"+probe.AreaError+";"+
//          coarseProbe.Background+";"+probe.Background+";"+
//          coarseProbe.Charge+";"+probe.Charge+";"+
//          coarseProbe.LowerMzBand+";"+probe.LowerMzBand+";"+
//          coarseProbe.Mz+";"+probe.Mz+";"+
//          coarseProbe.UpperMzBand+";"+probe.UpperMzBand+";"+
//          coarseProbe.LowerValley+";"+probe.LowerValley+";"+
//          coarseProbe.Peak+";"+probe.Peak+";"+
//          coarseProbe.UpperValley+";"+probe.UpperValley+";");
//      count++;
//      if (this.sameCgProbe(probe, coarseProbe))
//        inSameHash = true;
//    }
//    return inSameHash;
//  }
 
//  private boolean is3DProbeAlreadyUsed(Probe3D aProbe){
//    boolean alreadyUsed = false;
//    for (Probe3D probe: this.same3DProbes_){
//      if (this.sameCgProbe(probe, aProbe)&&
//          (probe.getEllipseMzPosition()-0.0001)<aProbe.getEllipseMzPosition()&&aProbe.getEllipseMzPosition()<(probe.getEllipseMzPosition()+0.0001)&&
//          (probe.getEllipseTimePosition()-0.001)<aProbe.getEllipseTimePosition()&&aProbe.getEllipseTimePosition()<(probe.getEllipseTimePosition()+0.001)&&
//          (probe.getEllipseMzStretch()-0.0001)<aProbe.getEllipseMzStretch()&&aProbe.getEllipseMzStretch()<(probe.getEllipseMzStretch()+0.0001)&&
//          (probe.getEllipseTimeStretch()-0.001)<aProbe.getEllipseTimeStretch()&&aProbe.getEllipseTimeStretch()<(probe.getEllipseTimeStretch()+0.001))
//        alreadyUsed = true;
//    }
//    return alreadyUsed;
//  }

  private boolean sameCgProbe(CgProbe probe1, CgProbe probe2)
  {
    boolean ok = false;
    /**** Here it is questionable if the CgProbe.Peak should be compared (for enhanced valley method) because the peak
     * can shift and then the area is painted green instead of red
     */
    if (probe1.Area >= (probe2.Area - 0.02)
        && probe1.Area <= (probe2.Area + 0.02)
        && probe1.AreaError >= (probe2.AreaError - 0.02)
        && probe1.AreaError <= (probe2.AreaError + 0.02)
        && probe1.Background >= (probe2.Background - 0.02)
        && probe1.Background <= (probe2.Background + 0.02)
        && probe1.Charge == probe2.Charge
        && probe1.LowerMzBand >= (probe2.LowerMzBand - 0.001)
        && probe1.LowerMzBand <= (probe2.LowerMzBand + 0.001)
        && probe1.LowerValley >= (probe2.LowerValley - 0.02)
        && probe1.LowerValley <= (probe2.LowerValley + 0.02)
        && probe1.Mz >= (probe2.Mz - 0.001)
        && probe1.Mz <= (probe2.Mz + 0.001)
        && probe1.Peak >= (probe2.Peak - 0.02)
        && probe1.Peak <= (probe2.Peak + 0.02)
        && probe1.UpperMzBand >= (probe2.UpperMzBand - 0.001)
        && probe1.UpperMzBand <= (probe2.UpperMzBand + 0.001)
        && probe1.UpperValley >= (probe2.UpperValley - 0.02)
        && probe1.UpperValley <= (probe2.UpperValley + 0.02)) {
      ok = true;
    }
    return ok;
  }
  
  private boolean sameProfileProbe(CgProbe probe1, CgProbe probe2)
  {
    boolean ok = false;
    /**** Here it is questionable if the CgProbe.Peak should be compared (for enhanced valley method) because the peak
     * can shift and then the area is painted green instead of red
     */
    if (probe1.Charge == probe2.Charge
        && probe1.LowerMzBand >= (probe2.LowerMzBand - 0.02)
        && probe1.LowerMzBand <= (probe2.LowerMzBand + 0.02)
        && probe1.LowerValley >= (probe2.LowerValley - 0.002)
        && probe1.LowerValley <= (probe2.LowerValley + 0.001)
        && probe1.Mz >= (probe2.Mz - 0.001)
        && probe1.Mz <= (probe2.Mz + 0.001)
        && probe1.Peak >= (probe2.Peak - 0.02)
        && probe1.Peak <= (probe2.Peak + 0.02)
        && probe1.UpperMzBand >= (probe2.UpperMzBand - 0.02)
        && probe1.UpperMzBand <= (probe2.UpperMzBand + 0.02)
        && probe1.UpperValley >= (probe2.UpperValley - 0.001)
        && probe1.UpperValley <= (probe2.UpperValley + 0.001)) {
      ok = true;
    }
    return ok;
  }
  
  /** clones the Vector before it unites the peaks*/
  private Vector<CgProbe> uniteTwinPeaksClone(@SuppressWarnings("rawtypes") Vector probes, LipidomicsChromatogram  chrom, int msLevel) throws CgException{
    Vector<Probe3D> clonedProbes = new Vector<Probe3D>();
    for (int i=0; i!=probes.size();i++){
      clonedProbes.add(new Probe3D((Probe3D)probes.get(i)));
    }
    return this.uniteTwinPeaks(clonedProbes, chrom, msLevel);
  }
  
  /** this is for the union of close peaks; this is done mainly because the other
   * isotopic peaks than the 0 isotope have no real change to be united and would
   * appear separated; but in general the procedure is quite similar to uniteTwinPeaks
   * @param probes
   * @param chrom
   * @return
   * @throws CgException
   */
  private Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> uniteClosePeaks(Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicAreas, Hashtable<Integer,LipidomicsChromatogram> chromHash, int msLevel)throws CgException{
    Hashtable<Integer,CgProbe> probes = new Hashtable<Integer,CgProbe>();
    for (Integer key : isotopicAreas.keySet()){
      Vector<CgProbe> areasOfOnePeak = isotopicAreas.get(key).get(0);
      if (areasOfOnePeak.size()>1) throw new CgException("The new edition with separated RTs has problems at the method uniteClosePeaks! Please report this case the developers");
      probes.put(key, areasOfOnePeak.get(0));
    } 
    if (probes.size()>0){
      Hashtable<Integer,Integer> excludeList = new Hashtable<Integer,Integer>();
      Vector<Integer> keys = new Vector<Integer>(probes.keySet());
      for (int i=0;i!=keys.size();i++){
        if (!excludeList.containsKey(i)){
          CgProbe probe1 = (CgProbe)probes.get(keys.get(i));
          for (int j=i+1;j<keys.size();j++){
            if (!excludeList.containsKey(j)){
              CgProbe probe2 = (CgProbe)probes.get(keys.get(j));
              // checks if the peaks are close to one another
              if (this.checkTimeSpace(probe1, probe2)<closePeakTimeTolerance_ && !isDeepValleyBetweenUnionPeaks(chromHash.get(0),probe1,probe2)){
                // Now the twin peaks are united; a new ellipse is calculated and intensities are extracted from the original file
                // a similar procedure like at detectPeak3D
                probe1 = uniteTwoPeaks(probe1,probe2,chromHash.get(0),excludeList,msLevel,j,null);
                // this unites peaks of other isotopes
                if (excludeList.containsKey(j)){
                  Hashtable<Integer,Vector<CgProbe>> isos1 = isotopicAreas.get(keys.get(i));
                  Vector<CgProbe> united = new Vector<CgProbe>();
                  united.add(probe1);
                  isos1.put(0, united);
                  Hashtable<Integer,Vector<CgProbe>> isos2 = isotopicAreas.get(keys.get(j));
                  int highestIsoNr = isos1.size();
                  if (isos2.size()>highestIsoNr) highestIsoNr = isos2.size();
                  for (int k=1; k<highestIsoNr; k++){
                    if (isos1.containsKey(k) && isos2.containsKey(k)){
                      if (isos1.get(k).size()>1) throw new CgException("The new edition with separated RTs has problems at the method uniteClosePeaks (higher isos)! Please report this case the developers");
                      if (isos2.get(k).size()>1) throw new CgException("The new edition with separated RTs has problems at the method uniteClosePeaks (higher isos)! Please report this case the developers");
                      Vector<CgProbe> united2 = new Vector<CgProbe>();
                      united2.add(uniteTwoPeaks(isos1.get(k).get(0),isos2.get(k).get(0),chromHash.get(k),new Hashtable<Integer,Integer>(),msLevel,j,probe1));
                      isos1.put(k, united2);
                    }else if (isos1.containsKey(k)){
                      isos1.put(k, isos1.get(k));
                    }else if (isos2.containsKey(k)){
                      isos1.put(k, isos2.get(k));
                    }
                  }
                  isotopicAreas.put(keys.get(i),isos1);
                  isotopicAreas.remove(keys.get(j)); 
                }
              }
            }
          }
        }
      }
    }
    return isotopicAreas;
  }
  
  private CgProbe uniteTwoPeaks (CgProbe probe1, CgProbe probe2,LipidomicsChromatogram chrom, Hashtable<Integer,Integer> excludeList, int msLevel,int element,CgProbe zeroIso) throws CgException{
    float oldLowerValley = probe1.LowerValley;
    float oldUpperValley = probe1.UpperValley;
    float oldPeak = probe1.Peak;
    float oldLowerMzBand = probe1.LowerMzBand;
    float oldUpperMzBand = probe1.UpperMzBand;

    boolean coversPeakTop = false;
    if (probe1.LowerValley<probe2.Peak&&probe2.Peak<probe1.UpperValley)coversPeakTop = true;
    if (probe2.LowerValley<probe1.Peak&&probe1.Peak<probe2.UpperValley)coversPeakTop = true;
    if (probe2.LowerValley<probe1.LowerValley)
      probe1.LowerValley = probe2.LowerValley;
    if (probe2.UpperValley>probe1.UpperValley)
      probe1.UpperValley = probe2.UpperValley;
    if (probe2.LowerMzBand>probe1.LowerMzBand)
      probe1.LowerMzBand = probe2.LowerMzBand;
    if (probe2.UpperMzBand>probe1.UpperMzBand)
      probe1.UpperMzBand = probe2.UpperMzBand;
    if (chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe2.Peak,chrom)][2]>chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe1.Peak,chrom)][2]){
      probe1.Peak = probe2.Peak;
    }
    float[] ellipse = Calculator.calculateEllipseFrom4Points(probe1.Peak, probe1.Mz-probe1.LowerMzBand,probe1.Peak, probe1.Mz+probe1.UpperMzBand, 
        probe1.LowerValley, probe1.Mz, probe1.UpperValley, probe1.Mz);
    if (ellipse[3]>profileMaxRange_)
      ellipse[3] = profileMaxRange_;
    CgProbe newProbe = probe1;
    if (probe1 instanceof Probe3D || probe2 instanceof Probe3D){
      CgProbe profileProbe = this.getProfileProbeWithoutQualityCheck(probe1,probe1.Mz-probe1.LowerMzBand,probe1.Mz+probe1.UpperMzBand,msLevel);
      Probe3D probe3D = new Probe3D(probe1,profileProbe,false,false,false,false,ellipse);
      LipidomicsChromatogram forQuant = new LipidomicsChromatogram(readJustIntensitiesOfInterest(probe3D,msLevel));
      forQuant.Background = probe1.Background;
      forQuant.getAreaRaw();
      probe3D.Area = forQuant.Area;
      newProbe = probe3D;
    }
//    System.out.println(probe3D.Mz+" Probe3D: "+probe3D.LowerValley+" ; "+probe3D.UpperValley);
    // this is for the special case if for whatever reason two times the same peak is returned
    if (((newProbe.Area*1.01)>probe1.Area&&(newProbe.Area*1.01)>probe2.Area)||coversPeakTop||zeroIso!=null){
      // this is a merging of 2 higher isotopic peaks, where probably a different peak is in between
      // then take the peak that is nearer to the +0 isotope
      if (zeroIso!=null && newProbe.Area>(2*(probe1.Area+probe2.Area))){
        float time1 = Math.abs(zeroIso.Peak-oldPeak);
        float time2 = Math.abs(zeroIso.Peak-probe2.Peak);
        if (time1>time2) probe1 = probe2;
        else{
          probe1.LowerValley = oldLowerValley;
          probe1.UpperValley = oldUpperValley;
          probe1.Peak = oldPeak;
          probe1.LowerMzBand = oldLowerMzBand;
          probe1.UpperMzBand = oldUpperMzBand;            
        }
      }else{
        probe1 = newProbe;
        excludeList.put(element,element);
      }
    }else{
      probe1.LowerValley = oldLowerValley;
      probe1.UpperValley = oldUpperValley;
      probe1.Peak = oldPeak;
      probe1.LowerMzBand = oldLowerMzBand;
      probe1.UpperMzBand = oldUpperMzBand;
    }
    return probe1;
  }
  

  /** sometimes ther is an artefact that peaks seam to be split (twin peaks)
   * this method unites this peaks
   * @param probes to be united
   * @param chrom the affected chromatogram
   * @return
   * @throws CgException
   */
  private Vector<CgProbe> uniteTwinPeaks(@SuppressWarnings("rawtypes") Vector probes, LipidomicsChromatogram chrom, int msLevel) throws CgException{
//  System.out.println("Uniting Twin peaks");
  Vector<CgProbe> unitedProbes = new Vector<CgProbe>();
//  Vector<CgProbe> localProbes = new Vector<CgProbe>(probes);
//  Vector<CgProbe> localUnitedProbes = new Vector<CgProbe>();
  if (probes.size()>0){
//    System.out.println("111111111111111111111");
//    while (localProbes.size()!=localUnitedProbes.size()){
//    if (localUnitedProbes!=null)
//      localProbes = new Vector<CgProbe>(localUnitedProbes);
//    localUnitedProbes = new Vector<CgProbe>();
    //excluded because already united to another Probe
    Hashtable<Integer,Integer> excludeList = new Hashtable<Integer,Integer>();
    // compare probes to one another
    for (int i=0;i!=probes.size();i++){
      if (!excludeList.containsKey(i)){
        Probe3D probe1 = (Probe3D)probes.get(i);
      for (int j=i+1;j<probes.size();j++){
        if (!excludeList.containsKey(j)){
          Probe3D probe2 = (Probe3D)probes.get(j);
//          System.out.println(probe1.LowerValley+";"+probe1.Peak+";"+probe1.UpperValley+";"+probe1.Mz+";"+probe1.AreaStatus);
//          System.out.println(probe2.LowerValley+";"+probe2.Peak+";"+probe2.UpperValley+";"+probe2.Mz+";"+probe2.AreaStatus);
//          System.out.println((probe1.getEllipseTimePosition()-probe1.getEllipseTimeStretch())+";"+(probe1.getEllipseTimePosition()+probe1.getEllipseTimeStretch()));
//          System.out.println((probe2.getEllipseTimePosition()-probe2.getEllipseTimeStretch())+";"+(probe2.getEllipseTimePosition()+probe2.getEllipseTimeStretch()));
          // checks if one probe is in the inner third around the peak-maximum of the other
//          System.out.println("2222222222222222222222222");
          boolean isInInnterThird = isOneProbeInOtherInnerThird(probe1,probe2,overlapPeakDistanceDivisor_,overlapFullDistanceDivisor_);
//          System.out.println("isInInnterThird: "+isInInnterThird);
          // checks if one probe is covered by the other one and if they are close by one another in 
          // m/z direction (since we have an additional degree of freedom)
          // OR if they are timely close to one another (isInInnerThird) it is checked if the ellipses overlap the Mz of the other peak
          if ((((probe2.LowerValley<=probe1.UpperValley&&probe1.UpperValley<=probe2.UpperValley)||
              (probe1.LowerValley<=probe2.UpperValley&&probe2.UpperValley<=probe1.UpperValley))&&
              ((probe1.Mz-0.002)<=probe2.Mz)&&(probe2.Mz<=(probe1.Mz+0.002)))||
              (isInInnterThird&&(probe1.getEllipseMzPosition()-probe1.getEllipseMzStretch())<=probe2.Mz)&&(probe2.Mz<=(probe1.getEllipseMzPosition()+probe1.getEllipseMzStretch()))){
//            System.out.println("111111111111");
//            System.out.println("333333333333333");
            if (!isDeepValleyBetweenPeaks(chrom,probe1,probe2)||isInInnterThird||
                (probe1.UpperValley>(probe2.Peak-(probe2.Peak-probe2.LowerValley)*0.66)&&probe1.LowerValley<probe2.Peak)||
                (probe2.UpperValley>(probe1.Peak-(probe1.Peak-probe1.LowerValley)*0.66)&&probe2.LowerValley<probe1.Peak)){
//              System.out.println("22222222222222");
              if (probe2.LowerValley<probe1.LowerValley)
                probe1.LowerValley = probe2.LowerValley;
              if (probe2.UpperValley>probe1.UpperValley)
                probe1.UpperValley = probe2.UpperValley;
              Probe3D useForProfilePoints = probe1;
              if (chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe2.Peak,chrom)][2]>chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe1.Peak,chrom)][2]){
                probe1.Peak = probe2.Peak;
                useForProfilePoints = probe2;
              }
              if (probe2.LowerMzBand>probe1.LowerMzBand)
                probe1.LowerMzBand = probe2.LowerMzBand;
              if (probe2.UpperMzBand>probe1.UpperMzBand)
                probe1.UpperMzBand = probe2.UpperMzBand;

              float[] ellipse = Calculator.calculateEllipseFrom4Points(probe1.Peak, probe1.Mz-probe1.LowerMzBand,
                probe1.Peak, probe1.Mz+probe1.UpperMzBand, 
                probe1.LowerValley, probe1.Mz, probe1.UpperValley, probe1.Mz);
              if (ellipse[3]>profileMaxRange_)
                ellipse[3] = profileMaxRange_;
              Probe3D probe3D = new Probe3D(probe1,useForProfilePoints,((Probe3D)probe1).isOverlapLowerMz(),
                ((Probe3D)probe1).isOverlapHigherMz(),((Probe3D)probe1).isOverlapBefore(),
                ((Probe3D)probe1).isOverlapLater(),ellipse);;
                LipidomicsChromatogram forQuant = new LipidomicsChromatogram(readJustIntensitiesOfInterest(probe3D,msLevel));
                forQuant.Background = probe1.Background; 
          
                forQuant.getAreaRaw();
                probe3D.Area = forQuant.Area;
                probe1 = probe3D;
                excludeList.put(j,j);
            }
          }
        }
      }
      unitedProbes.add(probe1); 
      }
    }
//    }  
  }
  return unitedProbes;
}
/*  
          if ((((probe2.LowerValley<=probe1.UpperValley&&probe1.UpperValley<=probe2.UpperValley)||
              (probe1.LowerValley<=probe2.UpperValley&&probe2.UpperValley<=probe1.UpperValley))&&
              ((probe1.Mz-twinPeakMzTolerance_)<=probe2.Mz)&&(probe2.Mz<=(probe1.Mz+twinPeakMzTolerance_)))||
              (isInInnerThird&&(probe1.getEllipseMzPosition()-probe1.getEllipseMzStretch())<=probe2.Mz)&&(probe2.Mz<=(probe1.getEllipseMzPosition()+probe1.getEllipseMzStretch()))){
//            System.out.println("111111111111");
            // if there is a deep valley between the peaks it can be assumed that they are not twin peaks
            // if there is a deep valley then at least the peaks have to be in the inner third to one another OR
            // the upper valley of the one peak has to be in the inner third or even more overlapping than the other
            if (!isDeepValleyBetweenPeaks(chrom,probe1,probe2)||isInInnerThird||
                (probe1.UpperValley>(probe2.Peak-(probe2.Peak-probe2.LowerValley)*0.66)&&probe1.LowerValley<probe2.Peak)||
                (probe2.UpperValley>(probe1.Peak-(probe1.Peak-probe1.LowerValley)*0.66)&&probe2.LowerValley<probe1.Peak)){
//              System.out.println("22222222222222");
              // Now the twin peaks are united; a new ellipse is calculated and intensities are extracted from the original file
              // a similar procedure like at detectPeak3D
              if (probe2.LowerValley<probe1.LowerValley)
                probe1.LowerValley = probe2.LowerValley;
              if (probe2.UpperValley>probe1.UpperValley)
                probe1.UpperValley = probe2.UpperValley;
              if (chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe2.Peak,chrom)][2]>chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe1.Peak,chrom)][2]){
                probe1.Peak = probe2.Peak;
              }
              if (probe2.LowerMzBand>probe1.LowerMzBand)
                probe1.LowerMzBand = probe2.LowerMzBand;
              if (probe2.UpperMzBand>probe1.UpperMzBand)
                probe1.UpperMzBand = probe2.UpperMzBand;

              float[] ellipse = Calculator.calculateEllipseFrom4Points(probe1.Peak, probe1.Mz-probe1.LowerMzBand,
                probe1.Peak, probe1.Mz+probe1.UpperMzBand, 
                probe1.LowerValley, probe1.Mz, probe1.UpperValley, probe1.Mz);
              Probe3D probe3D = new Probe3D(probe1,((Probe3D)probe1).isOverlapLowerMz(),
                ((Probe3D)probe1).isOverlapHigherMz(),((Probe3D)probe1).isOverlapBefore(),
                ((Probe3D)probe1).isOverlapLater(),ellipse);;
                LipidomicsChromatogram forQuant = this.readJustIntensitiesOfInterest(probe3D);
                forQuant.Background = probe1.Background; 
          
                forQuant.getAreaRaw();
                probe3D.Area = forQuant.Area;
                probe1 = probe3D;
                excludeList.put(j,j);
            }
          }
        }
      }
      unitedProbes.add(probe1); 
      }
    }
//    }  
  }
  
  
  return unitedProbes;
}*/

  /** This method checks if there is a deep valley between the two peaks
   * the range is defined by the the earlier and later peak  
   * for the calculation just smoothed values are used
   * first a highest and a lowest value are calculated
   * if the lowest value is higher than a threshold no deep-valley is assumed
   * @param chrom current chromatogram
   * @param probe1
   * @param probe2
   * @return isThereADeppValley
   */
  private boolean isDeepValleyBetweenPeaks(LipidomicsChromatogram chrom, CgProbe probe1, CgProbe probe2){
    boolean isDeepValley = true;
    float highestValue = 0f;
    float lowestValue = 10000000000f;
    int start = 0;
    int stop = 0;
    if (probe1.Peak<probe2.Peak){
      start = LipidomicsAnalyzer.findIndexByTime(probe1.Peak, chrom);
      stop = LipidomicsAnalyzer.findIndexByTime(probe2.Peak, chrom);
    }else{
      start = LipidomicsAnalyzer.findIndexByTime(probe2.Peak, chrom);
      stop = LipidomicsAnalyzer.findIndexByTime(probe1.Peak, chrom);
    }
    for (int i=start;i<stop;i++){
      if (chrom.Value[i][2]>highestValue)
        highestValue = chrom.Value[i][2];
      if (chrom.Value[i][2]<lowestValue)
        lowestValue = chrom.Value[i][2];
    }
    if (lowestValue>this.twinInBetweenCutoff_*highestValue)
      isDeepValley = false;
    return isDeepValley;
  }
  
  
  /** This method checks if there is a deep valley between the two peaks
   * the range is defined by the the earlier and later peak  
   * for the calculation just smoothed values are used
   * first a highest and a lowest value are calculated
   * if the lowest value is higher than a threshold no deep-valley is assumed
   * @param chrom current chromatogram
   * @param probe1
   * @param probe2
   * @return isThereADeppValley
   */
  private boolean isDeepValleyBetweenUnionPeaks(LipidomicsChromatogram chrom, CgProbe probe1, CgProbe probe2){
    boolean isDeepValley = true;
//    float highestValue = 0f;
    float lowestValue = Float.MAX_VALUE;
    int start = 0;
    int stop = 0;
    float highestProbe1 = 0;
    float highestProbe2 = 0;
    if (probe1.Peak<probe2.Peak){
      start = LipidomicsAnalyzer.findIndexByTime(probe1.Peak, chrom);
      stop = LipidomicsAnalyzer.findIndexByTime(probe2.Peak, chrom);
    }else{
      start = LipidomicsAnalyzer.findIndexByTime(probe2.Peak, chrom);
      stop = LipidomicsAnalyzer.findIndexByTime(probe1.Peak, chrom);
    }
    highestProbe1 = chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe1.Peak,chrom)][2];
    highestProbe2 = chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe2.Peak,chrom)][2];
    float highestValue = highestProbe1;
    if (highestProbe2<highestProbe1)
      highestValue = highestProbe2;
    for (int i=start;i<stop;i++){
      if (chrom.Value[i][2]<lowestValue)
        lowestValue = chrom.Value[i][2];
    }
    if (lowestValue>=unionInBetweenCutoff_*highestValue)
      isDeepValley = false;
    if (!isDeepValley) return isDeepValley;
    //this second check is to circumvent problems with the Savitzky-Golay smoothing
    lowestValue = Float.MAX_VALUE;
    highestProbe1 = chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe1.Peak,chrom)][1];
    highestProbe2 = chrom.Value[LipidomicsAnalyzer.findIndexByTime(probe2.Peak,chrom)][1];
    highestValue = highestProbe1;
    if (highestProbe2<highestProbe1)
      highestValue = highestProbe2;
    // this is when there not enough data points, the smoothed value may result on position where in the
    // raw data is no intensity
    if (highestValue==0){
       int startScan;
       int stopScan;
       if (probe1.Area<probe2.Area){
         if (probe1.Peak<probe2.Peak){
           startScan = start;
           if (probe1.UpperValley>probe2.Peak)
             stopScan = LipidomicsAnalyzer.findIndexByTime(probe1.Peak+(probe2.Peak-probe1.Peak)/2f, chrom);
           else
             stopScan = LipidomicsAnalyzer.findIndexByTime(probe1.UpperValley, chrom);
         }else{
           if (probe1.LowerValley<probe2.Peak)
             startScan = LipidomicsAnalyzer.findIndexByTime(probe1.Peak-(probe1.Peak-probe2.Peak)/2f, chrom);
           else
             startScan = LipidomicsAnalyzer.findIndexByTime(probe1.LowerValley, chrom);
           stopScan = stop;
         }
       }else{
         if (probe2.Peak<probe1.Peak){
           startScan = start;
           if (probe2.UpperValley>probe1.Peak)
             stopScan = LipidomicsAnalyzer.findIndexByTime(probe2.Peak+(probe1.Peak-probe2.Peak)/2f, chrom);
           else
             stopScan = LipidomicsAnalyzer.findIndexByTime(probe2.UpperValley, chrom);
         }else{
           if (probe2.LowerValley<probe1.Peak)
             startScan = LipidomicsAnalyzer.findIndexByTime(probe2.Peak-(probe2.Peak-probe1.Peak)/2f, chrom);
           else
             startScan = LipidomicsAnalyzer.findIndexByTime(probe2.LowerValley, chrom);
           stopScan = stop;
         }
       }
      for (int i=startScan;i<=stopScan;i++){
        if (chrom.Value[i][1]>highestValue)
          highestValue = chrom.Value[i][1];
      }      
    }
    
    for (int i=start;i<stop;i++){
      if (chrom.Value[i][1]<lowestValue)
        lowestValue = chrom.Value[i][1];
    }
    if (lowestValue>=unionInBetweenCutoff_*highestValue)
      isDeepValley = false;
    return isDeepValley;
  }
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  private Vector calculateIsoProbesWithinZeroIsotopeBoundaries(LipidomicsChromatogram chrom, Probe3D zeroProbe, Vector<CgProbe> otherAcceptedProbes, double multiplicationFactor, Vector<CgProbe> notAllowedIsoProbes,double zeroToOneMultiplicationFactor, boolean respectIntCutoff,int charge, int msLevel) throws CgException{
    CgProbe correctIsoProbe = null;
    Vector<Probe3D> allPossibleIsoProbes = new Vector<Probe3D>();
    Vector<CgProbe> notAllowedIsoProbe = new Vector<CgProbe>();
    int mainScan = LipidomicsAnalyzer.findIndexByTime(zeroProbe.Peak, chrom);
    // now a range where to look for the isotopic peak is defined by the peak borders of the previous chromatogram
    // TODO: here it would be better to really take the previous peak and not always the 0-isotope; this could prevent 
    // errors in +2 and higher chromatograms since it would narrow the range
    int lowestScan = LipidomicsAnalyzer.findIndexByTime(zeroProbe.LowerValley, chrom);
    int highestScan = LipidomicsAnalyzer.findIndexByTime(zeroProbe.UpperValley, chrom);
    // tolerance values around the threshold are calculated
    double lowerThreshold = multiplicationFactor*(1f/7f);
    double upperThreshold = multiplicationFactor*(2f/1f);
//    System.out.println("notAllowedIsoProbes: "+notAllowedIsoProbes.size());
    
    // the peak detection starts with the highest intensity of the main peak and 
    // goes then to the left and right side
    // the probes are addded to allPossibleIsoProbes
    CgProbe currentValidProbe = null;
    int currentScan = mainScan;
    Probe3D probe = null;
    try{
      probe = this.detectPeakThreeD(chrom, mainScan, respectIntCutoff,charge,msLevel);
      if (probe.AreaStatus == CgAreaStatus.OK){
        allPossibleIsoProbes.add(probe);
        currentValidProbe = probe;
        int lowerValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.LowerValley,chrom);
        if (lowerValleyIndex<currentScan)
          currentScan = lowerValleyIndex;
      }
    } catch (QuantificationException qcx){}

    currentScan--;
    while (currentScan>lowestScan){
      try{
        probe = this.detectPeakThreeD(chrom, currentScan, respectIntCutoff,charge,msLevel);
//      if (probe.AreaStatus == CgAreaStatus.AllreadyUsed)
//        probe.AreaStatus = CgAreaStatus.OK;
        if (probe.AreaStatus==CgAreaStatus.OK){
          this.checkDuplicate(allPossibleIsoProbes, probe);
          allPossibleIsoProbes.add(probe);
          currentValidProbe = probe;
        }
      } catch (QuantificationException qcx){}
      if (currentValidProbe!=null){
        int lowerValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.LowerValley,chrom);
        if (lowerValleyIndex<currentScan)
          currentScan = lowerValleyIndex;
      }
      currentScan--;
    }
    currentScan = mainScan;
    if (allPossibleIsoProbes.size()>0){
      int upperValleyIndex = LipidomicsAnalyzer.findIndexByTime(allPossibleIsoProbes.get(0).UpperValley,chrom);
        if (upperValleyIndex>currentScan)
          currentScan = upperValleyIndex;
    }
    currentScan++;
    while (currentScan<highestScan){
      try{
        probe = this.detectPeakThreeD(chrom, currentScan, respectIntCutoff,charge,msLevel);
//      if (probe.AreaStatus == CgAreaStatus.AllreadyUsed)
//        probe.AreaStatus = CgAreaStatus.OK;
        if (probe.AreaStatus==CgAreaStatus.OK){
          this.checkDuplicate(allPossibleIsoProbes, probe);
          allPossibleIsoProbes.add(probe);
          currentValidProbe = probe;
        }
      } catch (QuantificationException qcx){}
      if (currentValidProbe!=null){
        int upperValleyIndex = LipidomicsAnalyzer.findIndexByTime(currentValidProbe.UpperValley,chrom);
        if (upperValleyIndex>currentScan)
          currentScan = upperValleyIndex;
      }
      currentScan++;
    }
//    System.out.println("allPossibleIsoProbes.size(): "+allPossibleIsoProbes.size());
    // remove the duplicates
    Vector<Probe3D> removedDuplicates = new Vector<Probe3D>();
//    for (CgProbe notAllowedProbe : notAllowedIsoProbes){
//      System.out.println("notAllowed: "+notAllowedProbe.Area+";"+notAllowedProbe.LowerValley+";"+notAllowedProbe.Peak+";"+notAllowedProbe.UpperValley+";"+notAllowedProbe.AreaStatus);
//    }
    // now we check if the probe is possibly from a previous isotope
    for (Probe3D aProbe : allPossibleIsoProbes){
//      System.out.println(aProbe.Area+";"+aProbe.AreaStatus);
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        if (this.isFromPreviousIsotope(aProbe,notAllowedIsoProbes,zeroToOneMultiplicationFactor)){
          notAllowedIsoProbe.add(aProbe);
        }else
          removedDuplicates.add(aProbe);
      }  
    }
//    System.out.println("removedDuplicatesProbes.size(): "+removedDuplicates.size());
    // it could be that a peak has to be united to be corrected
    Vector<CgProbe> unitedPeaks = this.uniteTwinPeaksClone(removedDuplicates, chrom,msLevel);
    Vector<Probe3D> bothPeakTypes = new Vector<Probe3D>();
    for (Probe3D aProbe : removedDuplicates){
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        this.checkDuplicate(bothPeakTypes, aProbe);
        bothPeakTypes.add(aProbe);
      }
    }
    for (CgProbe aProbe : unitedPeaks){
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        this.checkDuplicate(bothPeakTypes, aProbe);
        bothPeakTypes.add((Probe3D)aProbe);
      }
    }
//    System.out.println("bothPeakTypes.size(): "+bothPeakTypes.size());
    // with the united peaks and the normal ones there are for sure duplicates since some are not united
    removedDuplicates = new Vector<Probe3D>();
    for (Probe3D aProbe : bothPeakTypes){
      if (aProbe.AreaStatus==CgAreaStatus.OK)
        removedDuplicates.add(aProbe);
    }
//    System.out.println("removedDuplicates: "+removedDuplicates.size());
    // now check which peak fulfills the area criterion
    Vector<Probe3D> fullfillsCriteria = new Vector<Probe3D>();
    float lowerArea = zeroProbe.Area*(float)lowerThreshold;
    float upperArea = zeroProbe.Area*(float)upperThreshold;
    for (Probe3D aProbe : removedDuplicates){
//      System.out.println(zeroProbe.Area+";"+multiplicationFactor+";"+lowerArea+";"+upperArea+";"+aProbe.Area);
      if (lowerArea<aProbe.Area&&aProbe.Area<upperArea)
        fullfillsCriteria.add(aProbe);
//      else if (lowerArea<aProbe.Area){
//        // this is just that analytes which are splitted are not per se excluded if they are splitted falsely
//        boolean foundAcceptableUnion = false;
//        for (CgProbe otherProbe : otherAcceptedProbes){
//          if (!otherProbe.isCoveredByThisProbe(zeroProbe)){
//            if ((otherProbe.LowerValley<zeroProbe.LowerValley&&otherProbe.UpperValley+0.01>=zeroProbe.LowerValley)||
//                (otherProbe.UpperValley>zeroProbe.UpperValley&&otherProbe.LowerValley-0.01<=zeroProbe.UpperValley)){
//              Probe3D newClonedProbe = new Probe3D(zeroProbe);
//              if (otherProbe.LowerValley<newClonedProbe.LowerValley)
//                newClonedProbe.LowerValley = otherProbe.LowerValley;
//              if (otherProbe.UpperValley>newClonedProbe.UpperValley)
//                newClonedProbe.UpperValley = otherProbe.UpperValley;
//              if (otherProbe.LowerMzBand>newClonedProbe.LowerMzBand)
//                newClonedProbe.LowerMzBand = otherProbe.LowerMzBand;
//              if (otherProbe.UpperMzBand>newClonedProbe.UpperMzBand)
//                newClonedProbe.UpperMzBand = otherProbe.UpperMzBand;
//              if (chrom.Value[LipidomicsAnalyzer.findIndexByTime(otherProbe.Peak,chrom)][2]>chrom.Value[LipidomicsAnalyzer.findIndexByTime(newClonedProbe.Peak,chrom)][2]){
//                newClonedProbe.Peak = otherProbe.Peak;
//              }
//              
//              float[] ellipse = Calculator.calculateEllipseFrom4Points(newClonedProbe.Peak, newClonedProbe.Mz-newClonedProbe.LowerMzBand,
//                  newClonedProbe.Peak, newClonedProbe.Mz+newClonedProbe.UpperMzBand, 
//                  newClonedProbe.LowerValley, newClonedProbe.Mz, newClonedProbe.UpperValley, newClonedProbe.Mz);
//              if (ellipse[3]>profileMaxRange_)
//                ellipse[3] = profileMaxRange_;
//              Probe3D probe3D = new Probe3D(newClonedProbe,false,false,false,false,ellipse);
//              LipidomicsChromatogram forQuant = this.readJustIntensitiesOfInterest(probe3D);
//              forQuant.Background = newClonedProbe.Background; 
//              forQuant.getAreaRaw();
//              probe3D.Area = forQuant.Area;
//              float newLowerArea = probe3D.Area*(float)lowerThreshold;
//              float newUpperArea = probe3D.Area*(float)upperThreshold;
////              System.out.println(aProbe.Area+";"+newLowerArea+";"+newUpperArea);
//              if (newLowerArea<aProbe.Area&&aProbe.Area<newUpperArea){
//                fullfillsCriteria.add(aProbe);
//                foundAcceptableUnion = true;
////                System.out.println("Now it seems to be OK!!");
//                break;
//              }
//            }
//          }
//        }
//        if (!foundAcceptableUnion)
//          notAllowedIsoProbe.add(aProbe);
      /*}*/else  
        notAllowedIsoProbe.add(aProbe);
    }
//    System.out.println("fullfillsCriteria: "+fullfillsCriteria.size());
    // if there are more than peaks that fulfill the area criterion we have to select which one to take
    if (fullfillsCriteria.size()>1){
      int position = this.findClosestToExpectedPosition(fullfillsCriteria,zeroProbe,multiplicationFactor);
      correctIsoProbe = fullfillsCriteria.get(position);
      fullfillsCriteria.remove(position);
      notAllowedIsoProbe.addAll(fullfillsCriteria);
    // if there is just one take it
    }else if (fullfillsCriteria.size()==1){
      correctIsoProbe = fullfillsCriteria.get(0);
    // if no one fulfills the criterion we try the 2D approach; if this still does not work we did not find any
    }else{
      CgProbe greedyProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
      // check if the greedy method works
      if (greedyProbe.AreaStatus==CgAreaStatus.OK && lowerArea<greedyProbe.Area&&greedyProbe.Area<upperArea){
        if (this.isFromPreviousIsotope(greedyProbe,notAllowedIsoProbes,zeroToOneMultiplicationFactor)){
          notAllowedIsoProbe.add(greedyProbe);
        }else{
          correctIsoProbe = greedyProbe;
        }
      }else{
        CgProbe standardProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.StandardValleyMethod,charge);
        if (standardProbe.Area<lowerArea)
          standardProbe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.EnhancedValleyMethod,charge);
//        System.out.println("At Standard Probe: "+lowerArea+";"+standardProbe.Area+";"+upperArea);
     // check if the standard valley method works
        if (lowerArea<standardProbe.Area&&standardProbe.Area<upperArea)
          if (!this.isFromPreviousIsotope(greedyProbe,notAllowedIsoProbes,zeroToOneMultiplicationFactor)){
            correctIsoProbe = standardProbe;
          }else{
            notAllowedIsoProbe.add(standardProbe);
          }
        else
          notAllowedIsoProbe.add(greedyProbe);
      }
    }
    Vector result = new Vector();
    result.add(correctIsoProbe);
    result.add(notAllowedIsoProbe);
//    System.out.println("returned not allowed: "+notAllowedIsoProbe.size());
    return result;
  }
  
  private int findClosestToExpectedPosition(Vector<Probe3D> probes, Probe3D zeroProbe, double multiplicationFactor){
    int correctProbePosition = -1;
//    Probe3D correctProbe = null;
    Hashtable<Integer,Integer> fromSameScoreToOriginalPosition = new Hashtable<Integer,Integer>();
    Vector<Integer> closeScore = new  Vector<Integer>();
    float peakDistance = zeroProbe.UpperValley-zeroProbe.LowerValley;
    int highestScore = 0;
    for (Probe3D probe: probes){
      int score = 0;
      if ((zeroProbe.Peak-peakDistance/6)<probe.Peak&&probe.Peak<(zeroProbe.Peak+peakDistance/6)){
//        System.out.println(1);
        score = score+3;
      }
      if ((zeroProbe.LowerValley-peakDistance/6)<probe.LowerValley&&probe.LowerValley<(zeroProbe.LowerValley+peakDistance/6)){
//        System.out.println(2);
        score = score+2;
      }else if (zeroProbe.LowerValley<probe.LowerValley&&probe.Peak<zeroProbe.UpperValley){
//        System.out.println(3);
        score = score+1;
      }
      if ((zeroProbe.UpperValley-peakDistance/6)<probe.UpperValley&&probe.UpperValley<(zeroProbe.UpperValley+peakDistance/6)){
//        System.out.println(4);
        score = score+2;
      }else if (zeroProbe.UpperValley>probe.UpperValley&&probe.Peak>zeroProbe.LowerValley ){
//        System.out.println(5);
        score = score+1;
      }
      closeScore.add(score);
      if (score>highestScore)
        highestScore = score;
//      System.out.println("score: "+score+" ; "+probe.Area);
    }
//    System.out.println("highestScore: "+highestScore);
    float idealArea = zeroProbe.Area*(float)multiplicationFactor;
    Vector<Probe3D> probesWithSameScore = new Vector<Probe3D>();
    float lowestAreaMultiplication = 100f;
    int lowestAreaPosition = 0;
    int currentAddPosition = 0;
    for (int i=0;i!=closeScore.size();i++){
      if (closeScore.get(i)==highestScore){
        Probe3D aProbe = probes.get(i);
        if (aProbe.Area>idealArea){
          if (aProbe.Area/idealArea<lowestAreaMultiplication){
            lowestAreaMultiplication = aProbe.Area/idealArea;
            lowestAreaPosition = currentAddPosition;
          }  
        }else{
          if (idealArea/aProbe.Area<lowestAreaMultiplication){
            lowestAreaMultiplication = idealArea/aProbe.Area;
            lowestAreaPosition = currentAddPosition;
          }
        }
        fromSameScoreToOriginalPosition.put(probesWithSameScore.size(),i);
        probesWithSameScore.add(probes.get(i));
        currentAddPosition++;
      }
    }
    if (probesWithSameScore.size()>1){
      correctProbePosition = fromSameScoreToOriginalPosition.get(lowestAreaPosition);
//      correctProbe = probesWithSameScore.get(lowestAreaPosition);
    }else{
      correctProbePosition = fromSameScoreToOriginalPosition.get(0);
//      correctProbe = probesWithSameScore.get(0);
    }
//    System.out.println("correctProbePosition: "+correctProbePosition);
    return correctProbePosition;
  }
  
  
  /** checks if the probe comes from another not-allowed isotope*/
  private boolean isFromPreviousIsotope(CgProbe aProbe,Vector<CgProbe> notAllowedIsoProbes,double multiplicationFactor){
    boolean isInNotAllowedIso = false;
    for (CgProbe otherIsoProbe: notAllowedIsoProbes){
      // 1.5 is the tolerance
      if ((otherIsoProbe.Area*multiplicationFactor*1.5)>aProbe.Area && isProbe1InProbe2Third(aProbe, otherIsoProbe,overlapPeakDistanceDivisor_,overlapFullDistanceDivisor_)){
        isInNotAllowedIso = true;
//        System.out.println("Is in not allowed1: "+otherIsoProbe.Area+";"+otherIsoProbe.LowerValley+";"+otherIsoProbe.Peak+";"+otherIsoProbe.UpperValley);
//        System.out.println("Is in not allowed2: "+aProbe.Area+";"+aProbe.LowerValley+";"+aProbe.Peak+";"+aProbe.UpperValley);
      }
        
    }
    return isInNotAllowedIso;
  }
  
  private static boolean isOneProbeInOtherInnerThird(CgProbe probe1, CgProbe probe2, float overlapPeakDistanceDivisor,
      float overlapFullDistanceDivisor){
    boolean isInThird = false;
    if (isProbe1InProbe2Third(probe1, probe2, overlapPeakDistanceDivisor, overlapFullDistanceDivisor))
      isInThird = true;
    if (isProbe1InProbe2Third(probe2, probe1, overlapPeakDistanceDivisor, overlapFullDistanceDivisor)) 
      isInThird = true;
    return isInThird;
  }
  
  private static boolean isProbe1InProbe2Third(CgProbe probe1, CgProbe probe2, float overlapPeakDistanceDivisor,
      float overlapFullDistanceDivisor){
    boolean isInInnerThird = false;
    if (probe1.Peak<probe2.Peak){
      float lowerThreshold1 = probe2.Peak-(probe2.Peak-probe2.LowerValley)/overlapPeakDistanceDivisor;
      float lowerThreshold2 = probe2.Peak-(probe2.UpperValley-probe2.LowerValley)/overlapFullDistanceDivisor;
//      System.out.println("lowerThreshold: "+lowerThreshold1+";"+lowerThreshold2+probe1.Peak);
      if (probe1.Peak>lowerThreshold1||probe1.Peak>lowerThreshold2){                    
        isInInnerThird = true;
      }  
    }else{
      float upperThreshold1 = probe2.Peak+(probe2.UpperValley-probe2.Peak)/overlapPeakDistanceDivisor;
      float upperThreshold2 = probe2.Peak+(probe2.UpperValley-probe2.LowerValley)/overlapFullDistanceDivisor;
//      System.out.println("upperThreshold: "+upperThreshold1+";"+upperThreshold2+";"+probe1.Peak);
      if (probe1.Peak<upperThreshold1||probe1.Peak<upperThreshold2){
        isInInnerThird = true;
      }  
    }
    return isInInnerThird;
  } 
  
  /** checks if the probe of a different isotopic pattern is near an accepted normal probe; 
   *  the ones which are far away are out of interest*/
  private boolean checkNearNormalProbe(CgProbe otherIsoProbe, Vector<CgProbe>probes,boolean[] isOfInterest){
    boolean nearNormalProbe = false;
    for (int i=0; i!= probes.size();i++){
      CgProbe probe = probes.get(i);
      if (this.checkTimeSpace(otherIsoProbe, probe)<isoNearNormalProbeTime_&&(isOfInterest==null || isOfInterest[i]))
        nearNormalProbe = true;
    }
    return nearNormalProbe;
  }
  
  private boolean otherIsoProbeInBetween(CgProbe highestProbe,CgProbe aProbe, Vector<CgProbe> otherIsoProbes){
    boolean isoInBetween = false;
    for (CgProbe otherIsoProbe:otherIsoProbes){
      if (otherIsoProbe.Area<(0.2f*aProbe.Area)) continue;
      if (highestProbe.Peak>aProbe.Peak){
        if (checkTimeSpace(highestProbe,aProbe)<isoInBetweenMaxTimeDistance_ && aProbe.Peak<otherIsoProbe.Peak&&otherIsoProbe.Peak<highestProbe.Peak)
          isoInBetween = true;
      }else{
        if (checkTimeSpace(highestProbe,aProbe)<isoInBetweenMaxTimeDistance_ && highestProbe.Peak<otherIsoProbe.Peak&&otherIsoProbe.Peak<aProbe.Peak)
          isoInBetween = true;      
      }
    }
    return isoInBetween;
  }
  
  private Vector<CgProbe> curateCurrentIsoProbes(Vector<CgProbe> currentIsoProbes, int isoNumber,Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicAreas,boolean[] isOfInterest){
    // curate the currentIsoProbes
    Vector<CgProbe> currentProbesOfInterest = new Vector<CgProbe>();
    for (Integer key : isotopicAreas.keySet()){
      if (isOfInterest[key]){
        if (isotopicAreas.get(key).containsKey(isoNumber)){
          for (CgProbe aProbe : isotopicAreas.get(key).get(isoNumber)){
            currentProbesOfInterest.add(aProbe);
          }
        }
      }
    }
    Vector<CgProbe> currentIsoProbesCurated = new Vector<CgProbe>();
    for (CgProbe probe1: currentIsoProbes){
      if (this.checkNearNormalProbe(probe1, currentProbesOfInterest, null)){
        boolean addIt = true;
        for (CgProbe probe2:currentProbesOfInterest){
          if (isOneProbeInOtherInnerThird(probe1, probe2,overlapPeakDistanceDivisor_,overlapFullDistanceDivisor_))
            addIt = false;
        }
        if (addIt)
          currentIsoProbesCurated.add(probe1);
      }
    }
    return currentIsoProbesCurated;
  }

  public void setGeneralBasePeakCutoff(float generalBasePeakCutoff)
  {
    this.generalBasePeakCutoff_ = generalBasePeakCutoff/1000f;
  }

  public float getHighestArea()
  {
    return highestArea_;
  }

  private void calculateNoiseCutoffValue(LipidomicsChromatogram chrom){
    Vector<Float> values = new Vector<Float>();
    for (int i=0; i!=chrom.Value.length;i++){
      values.add(chrom.Value[i][2]);
    }
    float median = Calculator.median(values);
    if (median>0){
      Vector<Float> doubleMedian = new Vector<Float>();
      for (int i =0; i!=values.size();i++){
        if (values.get(i)<(median*2)){
          doubleMedian.add(values.get(i));
        }
      }
      float[] floatArray = new float[doubleMedian.size()];
      for (int i=0; i!=doubleMedian.size();i++){
        floatArray[i] = doubleMedian.get(i);
      }
      float stdev = Calculator.stddeviation(floatArray);
      lowerIntensityThreshold_ = median;
      if (Float.isInfinite(stdev) || Float.isNaN(stdev))
        lowerIntensityThreshold_+=median;
      else
        lowerIntensityThreshold_+=noiseCutoffDeviationValue_*stdev;
    }
  }
  
  private void calculateMinIntThreshold(LipidomicsChromatogram chrom){
    float highestValue = 0f;
    for (int i=0; i!=chrom.Value.length;i++){
      if (chrom.Value[i][2]>highestValue) highestValue = chrom.Value[i][2];
    }
    if (lowerIntensityThreshold_<highestValue *this.minimumRelativeIntensity_) lowerIntensityThreshold_ = highestValue *this.minimumRelativeIntensity_.floatValue();
  }
  
  public CgProbe getProbeByManualSettings(LipidomicsChromatogram cr, boolean is3D, double startTime, double stopTime, double startMz, double stopMz, int charge, int msLevel)throws CgException{
    double timeCenter = (startTime+stopTime)/2d;
    cr.FindPeak(findIndexByTime((float)timeCenter, cr), 1, LipidomicsDefines.StandardValleyMethod);
    int lowValleyIndex =findIndexByTime((float)startTime, cr)+1;
    int upValleyIndex = findIndexByTime((float)stopTime, cr)+1;
    cr.LoValley = lowValleyIndex;
    cr.UpValley = upValleyIndex;
    if (cr.Peak<lowValleyIndex || cr.Peak>upValleyIndex)
      cr.Peak = (lowValleyIndex+upValleyIndex)/2;
    cr.GetBackground(50);
    cr.GetAreaAndTime();
    // this is a little bit a hack: the algorithm cannot integrate,if the peak
    // consists of one intensity, if this is the case, I set the background to zero
    // and integrate it (this is just important for FT data)
    if (cr.Area<=0){
      cr.Background = 0;
      cr.GetAreaAndTime();
    }
    //TODO: here this is just for charge 1
    CgProbe twoDProbe = new CgProbe(0,charge);
    twoDProbe = copyResults(twoDProbe, cr, false);
    CgProbe returnProbe = null;
    if (is3D){
      float broadChromSmoothRange = broadChromSmoothRange_*60f;
      double mz = (startMz+stopMz)/2d;
      double mzStretch = (stopMz-startMz)/2d;
      LipidomicsChromatogram broadChroma = new LipidomicsChromatogram(readAChromatogram((float)mz, (float)mzStretch, (float)mzStretch,msLevel,chromSmoothRange_,broadChromSmoothRepeats_,chromSmoothRange_,broadChromMeanSmoothRepeats_,(float)(timeCenter-broadChromSmoothRange),(float)(timeCenter+broadChromSmoothRange)));
      broadChroma.FindPeak(findIndexByTime((float)timeCenter,broadChroma), 1, LipidomicsDefines.GreedySteepnessPointImproved);
      broadChroma.LoValley = lowValleyIndex;
      broadChroma.UpValley = upValleyIndex;
      if (broadChroma.Peak<lowValleyIndex || broadChroma.Peak>upValleyIndex)
        broadChroma.Peak = (lowValleyIndex+upValleyIndex)/2;
      broadChroma.GetBackground(50);
      broadChroma.GetAreaAndTime();
      //TODO: here this is just for charge 1
      twoDProbe = new CgProbe(0,charge);
      twoDProbe = copyResults(twoDProbe, broadChroma, false);
      
      float[] ellipse = new float[4];
      ellipse[0] = (float)timeCenter;
      ellipse[1] = (float)mz;
      ellipse[2] = (float)((stopTime-startTime)/2d);
      ellipse[3] = (float)mzStretch;
      CgProbe profileProbe = getProfileProbeWithoutQualityCheck(twoDProbe, (float)startMz, (float)stopMz, msLevel);
      returnProbe = new Probe3D(twoDProbe,profileProbe,false,false,false,false,ellipse);
      LipidomicsChromatogram forQuant = new LipidomicsChromatogram(readJustIntensitiesOfInterest((Probe3D)returnProbe,msLevel));
      
      forQuant.Background = returnProbe.Background;    
      // just the raw intensities are extracted and stored as area; no smoothing is required
      forQuant.getAreaRaw();
      returnProbe.Area = forQuant.Area;
      if (returnProbe.Area>highestArea_)
        highestArea_ = returnProbe.Area;
      if (forQuant.getHighestIntensity()>highestIntensity_)
        highestIntensity_ = forQuant.getHighestIntensity();
    }else{
      returnProbe = twoDProbe;
    }
    returnProbe.AreaStatus = CgAreaStatus.OK;
    return returnProbe;
  }
  
  public Vector<Vector<CgProbe>> calculatePeakAtExactTimePosition(LipidParameterSet templateParam, /*int timeType,*/ int maxIsotope,int charge,int msLevel) throws CgException, QuantificationException {
    Vector<Vector<CgProbe>> results = new Vector<Vector<CgProbe>>();
    for (int i=0;i!=maxIsotope; i++){
      if (templateParam.getIsotopicProbes().size()>i){
        CgProbe templateProbe = templateParam.getIsotopicProbes().get(i).get(0);
        float mainRt = templateProbe.Peak;
        //if (timeType == LipidomicsDefines.MINUTES) mainRt = mainRt*60f;
        LipidomicsChromatogram chrom = new LipidomicsChromatogram(readAChromatogram(templateProbe.Mz, templateParam.LowerMzBand, templateParam.UpperMzBand, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        chrom.GetMaximumAndAverage();
        int mainScan = LipidomicsAnalyzer.findIndexByTime(mainRt, chrom);
        Vector<CgProbe> oneIso = new Vector<CgProbe>();
        CgProbe probe = null;
        try{probe = detectPeakThreeD(chrom,mainScan,false,charge,msLevel);}catch(QuantificationException qex){}
        if (probe!=null && probe.AreaStatus == CgAreaStatus.OK){
          oneIso.add(probe);
          results.add(oneIso);
          continue;
        }
        probe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.GreedySteepnessReductionMethod,charge);
        if (probe.AreaStatus == CgAreaStatus.OK){
          oneIso.add(probe);
          results.add(oneIso);
          continue;
        }
        probe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.EnhancedValleyMethod,charge);
        if (probe.AreaStatus == CgAreaStatus.OK){
          oneIso.add(probe);
          results.add(oneIso);
          continue;
        }
        probe = LipidomicsAnalyzer.calculateOneArea(chrom, mainScan, LipidomicsDefines.StandardValleyMethod,charge);
        if (probe.AreaStatus == CgAreaStatus.OK){
          oneIso.add(probe);
          results.add(oneIso);
          continue;
        }else
          break;
      }
    }
    return results;
  }
  
  public Vector<Vector<CgProbe>> calculatePeakAtExactProbePosition(LipidParameterSet templateParam, int maxIsotope,int charge, int msLevel) throws CgException {
    Vector<Vector<CgProbe>> results = new Vector<Vector<CgProbe>>();
    for (int i=0;i!=maxIsotope; i++){
      if (templateParam.getIsotopicProbes().size()>i){
        CgProbe probe = templateParam.getIsotopicProbes().get(i).get(0);
        LipidomicsChromatogram chrom = new LipidomicsChromatogram(readAChromatogram(probe.Mz, templateParam.LowerMzBand, templateParam.UpperMzBand, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        chrom.GetMaximumAndAverage();
        boolean is3D = false;
        float startMz = probe.Mz-probe.LowerMzBand;
        float stopMz = probe.Mz+probe.LowerMzBand;
        float startTime = probe.LowerValley;
        float stopTime = probe.UpperValley;
        if (probe instanceof Probe3D){
          is3D = true;
          Probe3D probe3D = (Probe3D)probe;
          startMz = probe3D.getEllipseMzPosition()-probe3D.getEllipseMzStretch();
          stopMz = probe3D.getEllipseMzPosition()+probe3D.getEllipseMzStretch();
          startTime = probe3D.getEllipseTimePosition()-probe3D.getEllipseTimeStretch();
          stopTime = probe3D.getEllipseTimePosition()+probe3D.getEllipseTimeStretch();
        }
        Vector<CgProbe> oneIso = new Vector<CgProbe>();
        CgProbe probeToAdd = getProbeByManualSettings(chrom, is3D, startTime, stopTime, startMz, stopMz,charge,msLevel);
        if (probeToAdd.Area>0){
          oneIso.add(probeToAdd);
          results.add(oneIso);
        }else break;
      }
    }
    return results;
  }
  
  public boolean justZeroValues(LipidomicsChromatogram profile){
    boolean justZero = true;
    for (int i=0; i!=profile.Value.length;i++){
      if (profile.Value[i][2]>profileIntThreshold_){
        justZero = false;
        break;
      }
    }   
    return justZero;
  }
  
  /**
   * 
   * @param px the CgProbe to enter the values
   * @param cx the chromatogram to get the values
   * @param ignoreZeros ignores empty values in the raw data in between - use this for profiles (not chromatograms)
   * @return
   */
  protected static CgProbe copyResults(CgProbe px, LipidomicsChromatogram cx, boolean ignoreZeros){
    px = Analyzer.copyResults(px, cx, ignoreZeros);
    px.greedyFragment = cx.isGreedyFragment();
    boolean loOverlap = false;
    if (cx.getLoSteepnessPoint()!=-1)
      loOverlap = true;
//    if (cx.Value[cx.LoValley][2]>(0.1*cx.Value[cx.LoValley][2]))
//      loOverlap = true;
    px.setLoOverlap(loOverlap);
    boolean upOverlap = false;
    if (cx.getUpSteepnessPoint()!=-1)
      upOverlap = true;
//    if (cx.Value[cx.UpValley][2]>(0.1*cx.Value[cx.UpValley][2]))
//      upOverlap = true;
    px.setUpOverlap(upOverlap);
//    px.loSteepnessPoint_ = cx.getLoSteepnessPoint();
//    px.upSteepnessPoint_ = cx.getUpSteepnessPoint();
    return px;
  }
  
  public static CgProbe calculateOneArea(CgChromatogram chrom,int scan, int mode, int charge){
    CgProbe probe = Analyzer.calculateOneArea(chrom, scan, mode);
    if (probe!=null) probe.Charge = charge;    
    return probe;
  }
  
  /**
   * calculates an area under the curve for MS2 fragments
   * @param fragmentMass the theoretical mass for the exctraction
   * @param fragmentFormula the chemical formula of the fragment
   * @param msLevel the level this fragment should be observed
   * @param charge the charge state of the fragment
   * @param isFromOtherSpecies defined if this fragment belongs to another isobaric species
   * @param probes the MS1 peak - defining the retention time limits
   * @return the calculated fragment area
   * @throws CgException thrown if there is something wrong with the chrom access
   */
  public CgProbe calculateMs2Area(double fragmentMass, String fragmentFormula, int msLevel, int charge, boolean isFromOtherSpecies, Vector<CgProbe> probes) throws CgException{
    //Juergen: the probes should belong to one peak, thus, only a range is extracted from the probes,
    // and all intensities inside the range are quantified 
    
    float[] startStopRt = this.getStartStopTimeFromProbes(probes);
    float lowestRt = startStopRt[0];
    float highestRt = startStopRt[1];

    //mzRange has to be read from the properties
    //this was the old code before the improved spectral caching
    //LipidomicsChromatogram chrom = new LipidomicsChromatogram(this.readAChromatogram((float)fragmentMass, msnMzTolerance_, msnMzTolerance_, msLevel,0f,0,0f,0,-1f,Float.MAX_VALUE));
    LipidomicsChromatogram chrom = new LipidomicsChromatogram(this.readJustIntensitiesOfInterest((float)fragmentMass-msnMzTolerance_,(float)fragmentMass+msnMzTolerance_,lowestRt,highestRt,msLevel));
    float area = 0f;
    CgProbe result = new CgProbe(0,charge,msLevel,fragmentFormula);
    float peakRt = 0f;
    float highestInt = 0f;
    result.AreaStatus = CgAreaStatus.NothingThere;
    if (chrom.getHighestIntensity()>0){
      int startScan = getStartPosition(chrom,lowestRt);
      for (int i=startScan; i!=chrom.Value.length; i++){
        float rt = chrom.Value[i][0];
        if (rt>highestRt) break;
        if (rt>lowestRt){
          area += chrom.Value[i][1];
          if (chrom.Value[i][1]>highestInt){
            highestInt = chrom.Value[i][1];
            peakRt = chrom.Value[i][0];
          }
        }
      }
    }
    if (area>0){
      result.AreaStatus = CgAreaStatus.OK;
      result.Area = area;
      result.AreaError = 0f;
      result.Background = 0f;
      result.Peak = peakRt;
      result.LowerValley = lowestRt;
      result.UpperValley = highestRt;
      result.Mz = (float)fragmentMass;
      result.LowerMzBand = msnMzTolerance_;
      result.UpperMzBand = msnMzTolerance_;
      result.isotopeNumber = 0;
      result.setFromOtherSpecies(isFromOtherSpecies);
    }
//    for (CgProbe probe : probes){
//      System.out.println(probe.LowerValley+";"+probe.UpperValley);
//    }
//    String basePath = "E:\\lipidomicsMS2\\20130729\\";
//    printChromaToFile(chrom,basePath+String.valueOf(fragmentMass)+".png",1);
    return result;
  }
  
  /**
   * returns the retention times of MSn spectra that are within the detected MS1 peaks
   * @param msLevel the MS-level of the spectra
   * @param probes the detected MS1 peaks
   * @return the sorted retention times of MSn spectra that are within the detected MS1 peaks; first key scan number; value: retention time
   */
  public LinkedHashMap<Integer,Float> getMSnSpectraRetentionTimes(int msLevel, Vector<CgProbe> probes){
    LinkedHashMap<Integer,Float> rts = new LinkedHashMap<Integer,Float>();
    float[] startStopRt = this.getStartStopTimeFromProbes(probes);
    float lowestRt = startStopRt[0];
    float highestRt = startStopRt[1];
    Hashtable<Integer,Hashtable<Integer,String>> spectraCache = getMSnSpectraCache();
    if (!spectraCache.containsKey(msLevel)) return rts;
    Hashtable<Integer,String> spectra = spectraCache.get(msLevel);
    Hashtable<Integer,Float> retTimes = reader_.getRetentionTimes(msLevel);
    List<Integer> consScanNumbers = new ArrayList<Integer>(spectra.keySet());
    Collections.sort(consScanNumbers);
    for (Integer consScanNumber : consScanNumbers){
      float rt = retTimes.get(consScanNumber);
      int scanNumber = Integer.parseInt(spectra.get(consScanNumber).substring(0,spectra.get(consScanNumber).indexOf(" ")));
      if (lowestRt<rt && rt<=highestRt)
        rts.put(scanNumber, rt);
    }
    return rts;
  }
  
  private int getStartPosition(CgChromatogram chrom, float startRt){
    int startScan = chrom.Value.length;
    for (int j=0; j<chrom.Value.length; j+=100){
      if (j+100>=chrom.Value.length || chrom.Value[j+100][0]>startRt){
        for (int k=j; k!=j+100; k+=10){
          if (k+10>=chrom.Value.length || chrom.Value[k+10][0]>startRt){
            for (int i=k; i<chrom.Value.length; i++)
            {
              if (i+1>=chrom.Value.length || chrom.Value[i+1][0]>startRt) return i;
            }
          }  
        }
      }
    }
    return startScan;

  }
  
  public Hashtable<Integer,Float> extractBasePeakValues(Vector<Integer> levels, Vector<CgProbe> probes) throws CgException {
    return  extractBasePeakValues(levels, probes, msnMzTolerance_);  
  }
  
  /**
   * calculates extraction ranges for found scans
   * @param spectrLevelsRequired which spectral ranges are required for this set of rules
   * @return start stop times in the form of a "Range"
   */
  public Vector<Range> findSingleSpectraRanges(int[] spectrLevelsRequired){
    Vector<Range> ranges = new Vector<Range>();
    Hashtable<Integer,Hashtable<Integer,String>> spectraCache = getMSnSpectraCache();
    Hashtable<Integer,List<Float>> rtsSorted = new Hashtable<Integer,List<Float>>();
    Vector<Integer> foundLevels = new Vector<Integer>();
    for (int i=spectrLevelsRequired[0]; i<(spectrLevelsRequired[1]+1); i++){
      List<Float> rts = new ArrayList<Float>();
      if (spectraCache.containsKey(i)){
        Hashtable<Integer,String> spectra = spectraCache.get(i);
        Hashtable<Integer,Float> retTimes = getRetentionTimes(i);
        for (Integer consScanNumber : spectra.keySet()) rts.add(retTimes.get(consScanNumber));
        Collections.sort(rts);
        rtsSorted.put(i, rts);
        foundLevels.add(i);
      }
    }
    if (foundLevels.size()==0) return ranges; 
    float startTime = 0f;
    float stopTime = 0f;
    List<Float> rts = rtsSorted.get(foundLevels.get(0));
    Hashtable<Integer,Integer> levelStopped = new Hashtable<Integer,Integer>();
    for (int level : foundLevels) levelStopped.put(level, 0);
    for (int i=0; i<rts.size(); i++){
      startTime = stopTime;
      if (i==0) startTime = rts.get(i)-rts.get(i)/10000f;
      float firstValue = rts.get(i);
      float secondValue = 0f;
      if ((i+1)<rts.size()) secondValue = rts.get(i+1);
      else secondValue = Float.MAX_VALUE;
      for (int j=1; j<foundLevels.size(); j++){
        int level = foundLevels.get(j);
        List<Float> lowRts = rtsSorted.get(level);
        int stopped = levelStopped.get(level);
        for (int k=stopped; k<lowRts.size(); k++){
          float rt = lowRts.get(k);
          stopped = k;
          if (firstValue<rt && rt<secondValue) firstValue = rt;
          if (rt>secondValue) break;
        }
        levelStopped.put(level, stopped);
      }
      if ((i+1)==rts.size()) secondValue = firstValue*1.001f;
      stopTime = (firstValue+secondValue)/2f;
      ranges.add(new Range(startTime,stopTime));
    }
    return ranges;
  }

  /**
   * merges several MS identifcations into a single one
   * @param hits the MS identifications to be merged
   * @return the single identificatin for each isotope
   * @throws CgException if there is something wrong with the extraction
   */
  public Hashtable<Integer,Vector<CgProbe>> mergeIdentifications(Collection<LipidParameterSet> hits) throws CgException{
    Hashtable<Integer,Vector<Float>> maxDimensions = extractMaxDimensions(hits);
    boolean negative = false;
    if (maxDimensions.size()>0 && maxDimensions.containsKey(-1)) negative = true;
    Hashtable<Integer,Vector<CgProbe>> results = new Hashtable<Integer,Vector<CgProbe>>();
    for (int i=0; i!=maxDimensions.size(); i++){
      Vector<Float> dimensions = null;
      int isotopeNumber = i;
      if (negative) isotopeNumber = isotopeNumber*-1;
      dimensions = maxDimensions.get(isotopeNumber);
      float lowestMz = dimensions.get(0);
      float highestMz = dimensions.get(1);
      float lowestRt = dimensions.get(2);
      float highestRt = dimensions.get(3);
      float[] ellipse = new float[]{(lowestRt+highestRt)/2f,(lowestMz+highestMz)/2f,(highestRt-lowestRt)/2f,(highestMz-lowestMz)/2f};
      CgProbe probe = new CgProbe(-1,hits.iterator().next().getCharge());
      LipidomicsChromatogram chrom = new LipidomicsChromatogram(readAChromatogram(ellipse[1], ellipse[3]/2, ellipse[3]/2, 1, chromSmoothRange_,chromSmoothRepeats_));
      chrom.GetMaximumAndAverage();
      chrom.detectPeakIntensityInsideBorders(LipidomicsAnalyzer.findIndexByTime(lowestRt,chrom),LipidomicsAnalyzer.findIndexByTime(highestRt,chrom));
      chrom.GetBackground(50);
      chrom.GetAreaAndTime();
      Analyzer.copyResults(probe, chrom, false);
      probe.AreaStatus = CgAreaStatus.OK;
      CgProbe profileProbe = getProfileProbeWithoutQualityCheck(probe, lowestMz, highestMz, probe.getMsLevel());
      Probe3D probe3D = new Probe3D(probe,profileProbe,false,false,false,false,ellipse);
      
      LipidomicsChromatogram forQuant = new LipidomicsChromatogram(readJustIntensitiesOfInterest(probe3D,1));

      
      forQuant.Background = probe.Background; 
      forQuant.getAreaRaw();
      probe3D.Area = forQuant.Area;
      Vector<CgProbe> probes = new Vector<CgProbe>();
      probes.add(probe3D);
      results.put(i, probes);
    }
    return results;
  }
  
  /**
   * extracts the minimum and maximum m/z and RT values of several identifications
   * @param hits the MS1 identifications
   * @return a hashtable containing the minimum and maximum m/z and RT values;
   * the key: the isotope number; value: max values - 1. lowestMz; 2. highestMz; 3. lowest RT; 4. highest RT
   * 
   */
  private Hashtable<Integer,Vector<Float>> extractMaxDimensions(Collection<LipidParameterSet> hits){
    Hashtable<Integer,Vector<Float>> maxDimensions = new Hashtable<Integer,Vector<Float>>();
    int maxIsotope = 0;
    for (LipidParameterSet set : hits){
      if (set.getIsotopicProbes().size()>maxIsotope) maxIsotope = set.getIsotopicProbes().size();
    }
    for (int i=0; i!=maxIsotope; i++){
      float lowestMz = Float.MAX_VALUE;
      float highestMz = 0f;
      float lowestRt = Float.MAX_VALUE;
      float highestRt = 0f;
      int isotopeNumber = 0;
      for (LipidParameterSet set : hits){
        if ((set.getIsotopicProbes().size()-1)<i) continue;
        Vector<CgProbe> probes = set.getIsotopicProbes().get(i);
        for (CgProbe probe : probes){
          isotopeNumber = probe.isotopeNumber;
          if ((probe.Mz-probe.LowerMzBand)<lowestMz) lowestMz = probe.Mz-probe.LowerMzBand;
          if ((probe.Mz+probe.UpperMzBand)>highestMz) highestMz = probe.Mz+probe.UpperMzBand;
          if (probe.LowerValley<lowestRt) lowestRt = probe.LowerValley;
          if (probe.UpperValley>highestRt) highestRt = probe.UpperValley;
        }
      }
      Vector<Float> dimensions = new Vector<Float>();
      dimensions.add(lowestMz);
      dimensions.add(highestMz);
      dimensions.add(lowestRt);
      dimensions.add(highestRt);
      maxDimensions.put(isotopeNumber, dimensions);
    }
    return maxDimensions;
  }

  /**
   * 
   * @return the area cutoff relative to the highest peak identified for this m/z
   */
  public float getRelativeAreaCutoff()
  {
    return relativeAreaCutoff_;
  }

  /**
   * 
   * @return the area cutoff of peaks farer away from the highest peak identified for this m/z
   */
  public float getRelativeFarAreaCutoff()
  {
    return relativeFarAreaCutoff_;
  }
  
  
  /**
   * 
   * @return area multiplication factor to remove very small uninteresting peaks
   */
  public int getPeakDiscardingAreaFactor()
  {
    return peakDiscardingAreaFactor_;
  }

  /**
   * sets to relative area cutoff for peaks of one chromatogram
   * @param relativeAreaCutoff area cutoff relative to the highest peak identified for this m/z
   * @param relativeFarAreaCutoff area cutoff of peaks farer away from the highest peak identified for this m/z
   * @param peakDiscardingAreaFactor area multiplication factor to remove very small uninteresting peaks
   * @param overWriteDiscardingFactorOnlyIfBigger if true, the discarding factor will only be overwritten if it is bigger
   *                                              (alllowing smaller values); if false, it will be overwritten anyway
   */
  public void setAreaCutoffs(float relativeAreaCutoff, float relativeFarAreaCutoff,
      int peakDiscardingAreaFactor, boolean overWriteDiscardingFactorOnlyIfBigger){
    this.relativeAreaCutoff_ = relativeAreaCutoff;
    this.relativeFarAreaCutoff_ = relativeFarAreaCutoff;
    if (!overWriteDiscardingFactorOnlyIfBigger || peakDiscardingAreaFactor>this.peakDiscardingAreaFactor_)
      this.peakDiscardingAreaFactor_ = peakDiscardingAreaFactor; 
  }
  
  
  /**
   * this method checks if any of the identified 0 isotopic peaks share the same peak center
   * @param set1 the first MS1 identification
   * @param set2 the second MS1 identification
   * @return true if there any of the peak centers the same
   */
  public static boolean isPeakCenterTheSame(LipidParameterSet set1, LipidParameterSet set2){
    boolean theSame = false;
    for (CgProbe probe1 : set1.getIsotopicProbes().get(0)){
      for (CgProbe probe2 : set2.getIsotopicProbes().get(0)){
        if (isOneProbeInOtherInnerThird(probe1,probe2,LipidomicsConstants.getOverlapPeakDistanceDivisor(),LipidomicsConstants.getOverlapFullDistanceDivisor()))
          return true;
      }
    }
    return theSame;
  }
  
  /**
   * returns one chromatogram at a speciefied m/z value and MS level
   * @param mz the m/z value for the chromatogram
   * @param msLevel the MS level
   * @return the chromatogram
   * @throws CgException thrown if something is wrong in the reading procedure form the chrom file
   */
  public LipidomicsChromatogram readOneChromatogram(float mz, int msLevel) throws CgException{
    return new LipidomicsChromatogram(readAChromatogram(mz, coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
  }
  
  /**
   * returns peak borders, where the intensities of the peak drop below a certain percentual threshold
   * @param percent the percentage for the threshold
   * @param chrom the chromatogram
   * @param sets the identifications where the new percentual threshold shall be calculate
   * @return float[0] = border start value; float[1] = border stop value
   */
  public float[] calculatePercentualBorders(float percent, LipidomicsChromatogram chrom, Vector<LipidParameterSet> sets){
    float start = Float.MAX_VALUE;
    float stop = 0f;
    CgProbe intenseProbe = null;
    for (LipidParameterSet set : sets){
      
      for (CgProbe probe: set.getIsotopicProbes().get(0)){
        if (intenseProbe==null || probe.Area>intenseProbe.Area) intenseProbe = probe;
      }
    }
    if (intenseProbe!=null){
      float[] borders = chrom.findPercentualBorders(percent,intenseProbe.LowerValley,intenseProbe.UpperValley);
      start = borders[0];
      stop = borders[1];
    }
    
    float[] result = new float[2];
    result[0] = start;
    result[1] = stop;
    return result;
  }

  public float getMsnMzTolerance()
  {
    return msnMzTolerance_;
  }
  
  public boolean getUseCuda()
  {
	return useCuda_;
  }
  
  public SavGolJNI getSavGolJNI()
  {
	  return sav_gol_jni_;
  }
  
  /**
   * in the isobaric peak separation procedure, hard limits for the peak are set
   * this method calculates the peak areas only inside these hard limits
   * @param set the set holding the information about the previous identification
   * @param msLevel the MS level of the identification (usually MS1)
   * @param peak the newly calculated peak center (out of MSn)
   * @throws CgException if there is something wrong with the chrom access
   */
  public void recalculatePeaksAccordingToHardLimits(LipidParameterSet set, int msLevel, float peak) throws CgException{
    Vector<Vector<CgProbe>> newProbes = new Vector<Vector<CgProbe>>();
    float totalArea = 0f;
    for (int i=0; i!=set.getIsotopicProbes().size(); i++){
      Vector<CgProbe> probes = set.getIsotopicProbes().get(i);
      Vector<CgProbe> newProbesOfIso = new Vector<CgProbe>();
      LipidomicsChromatogram chrom = null;
      if (chromHash_.containsKey(i)) chrom = chromHash_.get(i);
      else{
        chrom = new LipidomicsChromatogram(readAChromatogram(set.Mz[0], coarseChromMzTolerance_, coarseChromMzTolerance_, msLevel, chromSmoothRange_,chromSmoothRepeats_));
        chromHash_.put(i, chrom);
      }
      for (CgProbe aProbe: probes){
        float lowerLimit = set.getLowerRtHardLimit();
        float upperLimit = set.getUpperRtHardLimit();
        if (lowerLimit>=0){
          if (lowerLimit>=aProbe.UpperValley) continue; 
        } else if (upperLimit>=0){
          if (upperLimit<=aProbe.LowerValley) continue;
        }
        CgProbe newProbe = aProbe;
        if ((lowerLimit>=0 && newProbe.LowerValley<lowerLimit) ||
            (upperLimit>=0 && newProbe.UpperValley>upperLimit)){
          float loVal = newProbe.LowerValley;
          float upVal = newProbe.UpperValley;
          if (lowerLimit>=0 && newProbe.LowerValley<lowerLimit)
            loVal = lowerLimit;
          if (upperLimit>=0 && newProbe.UpperValley>upperLimit)
            upVal = upperLimit;
          if (newProbe instanceof Probe3D){
            Probe3D probe3D = new Probe3D((Probe3D)newProbe);
            if (lowerLimit>=0 && newProbe.LowerValley<lowerLimit)
              probe3D.setLowerHardRtLimit(lowerLimit);
            if (upperLimit>=0 && newProbe.UpperValley>upperLimit)
              probe3D.setUpperHardRtLimit(upperLimit);
            LipidomicsChromatogram forQuant = new LipidomicsChromatogram(readJustIntensitiesOfInterest(probe3D,msLevel));
            forQuant.Background = probe3D.Background;    
            // just the raw intensities are extracted and stored as area; no smoothing is required
            forQuant.getAreaRaw();
            probe3D.Area = forQuant.Area;
            chrom.detectPeakIntensityInsideBorders(LipidomicsAnalyzer.findIndexByTime(loVal,chrom), LipidomicsAnalyzer.findIndexByTime(upVal,chrom));
            probe3D.LowerValley = loVal;
            probe3D.UpperValley = upVal;
            newProbe = probe3D;
          }else{
            newProbe = new CgProbe(newProbe);
            chrom.detectPeakIntensityInsideBorders(LipidomicsAnalyzer.findIndexByTime(loVal,chrom), LipidomicsAnalyzer.findIndexByTime(upVal,chrom));
            chrom.GetBackground(50);
            chrom.GetAreaAndTime();
            Analyzer.copyResults(newProbe, chrom, false);
            newProbe.AreaStatus = CgAreaStatus.OK;
          }
          newProbe.Peak = peak;
          if (newProbe.Peak<loVal) newProbe.Peak = loVal;
          if (newProbe.Peak>upVal) newProbe.Peak = upVal;
          newProbesOfIso.add(newProbe);
        }else {
          newProbesOfIso.add(newProbe);
        }
      }
      if (newProbesOfIso.size()>0){
        newProbes.add(newProbesOfIso);
        for (CgProbe probe : newProbesOfIso) totalArea += probe.Area;
      } else {
        break;
      }
    }
    if (newProbes.size()==0)
      throw new CgException("A peak split cannot be performed because one of the partners does not deliver MS1 peak areas!");
    set.setIsotopicProbes(newProbes);
    set.Area = totalArea;
  }
  
  private void checkIfSmoothingDeviatedPeakPosition(LipidomicsChromatogram chrom, CgProbe probe){
    int peak = findIndexByTime(probe.Peak, chrom);
    if (chrom.Value[peak][1]>0) return;
    int start = findIndexByTime(probe.LowerValley, chrom);
    int stop = findIndexByTime(probe.UpperValley, chrom);
    float highestInt = 0f;
    int newPeak = -1;
    for (int i=start; i<=stop; i++){
      if (chrom.Value[i][1]>highestInt){
        highestInt = chrom.Value[i][1];
        newPeak = i;
      }
    }
    if (newPeak>-1) probe.Peak = chrom.Value[newPeak][0];
//    System.out.println("!!!!!!!!!!!!!!!!!!!!!! "+chrom.Value[peak][1]+";"+chrom.Value[start][0]/60f+";"+chrom.Value[stop][0]/60f);
//    int lowerCount = peak;
//    int upperCount = peak;
//    boolean search = true;
//    while (search){
//      boolean takeUpperCount = false;
//      boolean takeLowerCount = false;
//      if ((lowerCount-1)<start||lowerCount==0) takeUpperCount=true;
//      if (((upperCount+1)>stop || (upperCount+1)==chrom.Value.length)) takeLowerCount=true;
//      if (takeUpperCount && takeLowerCount) break;
//      if (!takeUpperCount && !takeLowerCount){
//        if (peak-lowerCount==upperCount-peak){
//          System.out.println(chrom.Value[lowerCount-1][1]+";"+chrom.Value[upperCount+1][1]);
//          if (chrom.Value[lowerCount-1][1]>chrom.Value[upperCount+1][1]) takeLowerCount = true;
//          else takeUpperCount = true;
//        } else if (peak-lowerCount<upperCount-peak){
//          takeLowerCount=true;
//        }else{
//          takeUpperCount = true;
//        }
//      }
//      if (takeLowerCount){
//        lowerCount--;
//        if (chrom.Value[lowerCount][1]>0){
//          probe.Peak = chrom.Value[lowerCount][0];
//          search = false;
//        }
//      }
//      if (takeUpperCount){
//        upperCount++;
//        if (chrom.Value[upperCount][1]>0){
//          probe.Peak = chrom.Value[upperCount][0];
//          search = false;
//        }
//      }
//    }
//    System.out.println(lowerCount+";"+upperCount+";"+peak);
  }
  
  /**
   * detects an m/z profile for a peak that was found by chromatogram analysis
   * @param aProbe the detected peak in a chromatogram
   * @param msLevel the MS-level of the extracted profile
   * @param smooth should the profile be smoothed
   * @return an m/z profile for a peak that was found by chromatogram analysis
   * @throws CgException if there is something wrong with the chrom access
   */
  private LipidomicsChromatogram readASingleProfile(CgProbe aProbe, int msLevel, boolean smooth) throws CgException{
    Vector<CgProbe> theProbes = new Vector<CgProbe>();
    theProbes.add(aProbe);
    // here the raw profile is extracted from the file and smoothed 
    float mzRange = this.profileMzRange_;
    Vector<CgChromatogram> profiles = this.readProfiles(theProbes, mzRange, this.profileTimeTolerance_,5f,
      profileSmoothRange_,profileSmoothRepeats_, profileSmoothRange_, profileMeanSmoothRepeats_, msLevel, smooth);
    LipidomicsChromatogram profile = new LipidomicsChromatogram(profiles.get(0));
    if (justZeroValues(profile)){
      profiles = this.readProfiles(theProbes, mzRange, this.broaderProfileTimeTolerance_,5f,
          profileSmoothRange_,profileSmoothRepeats_, profileSmoothRange_, profileMeanSmoothRepeats_, msLevel, smooth);
      profile = new LipidomicsChromatogram(profiles.get(0));
    }
    return profile;
  }
  
  /**
   * detects an m/z profile peak without checking if the peak is OK - used for manual setting of the peak borders for getting 10% and 50% values relative to the apex
   * @param probe the MS1-peak containing the borders
   * @param startMz the start m/z value
   * @param stopMz the stop m/z value
   * @param msLevel the MS-level of the extracted profile
   * @return an m/z profile peak without checking if the peak is OK - used for manual setting of the peak borders for getting 10% and 50% values relative to the apex
   * @throws CgException
   */
  private CgProbe getProfileProbeWithoutQualityCheck(CgProbe probe, float startMz, float stopMz, int msLevel) throws CgException{
    LipidomicsChromatogram profile = readASingleProfile(probe,msLevel,false);
    profile.LoValley = LipidomicsAnalyzer.findIndexByTime((float)startMz,profile);
    profile.UpValley = LipidomicsAnalyzer.findIndexByTime((float)stopMz,profile);
    profile.Good = true;
    profile.detectPeakIntensityInsideBorders(LipidomicsAnalyzer.findIndexByTime((float)startMz,profile),LipidomicsAnalyzer.findIndexByTime((float)stopMz,profile));
    CgProbe profileProbe = new CgProbe(0,probe.Charge);
    profileProbe = LipidomicsAnalyzer.copyResults(profileProbe, profile, true);
    return profileProbe;
  }
  
  /**
   * caches the spectra of a certain precursor mass
   * @param from lower m/z threshold for precursor m/z
   * @param to upper m/z threshold for precursor m/z
   * @param from lower retention time
   * @param top upper retention time
   * @return
   * @throws CgException
   */
  public Hashtable<Integer,Boolean> prepareMSnSpectraCache(float startMz, float stopMz, float startTime, float stopTime) throws CgException{
    return reader_.prepareMSnSpectraCache(startMz, stopMz, startTime, stopTime,LipidomicsConstants.getMs2MinIntsForNoiseRemoval());
  }

  /**
   * checks in the MSn spectra cache which MSn levels are present
   * @return
   */
  public Hashtable<Integer,Boolean> checkMSnLevels(){
    Hashtable<Integer,Hashtable<Integer,String>> cache = reader_.getMSnSpectraCache();
    Hashtable<Integer,Boolean> availableLevels = new Hashtable<Integer,Boolean>();
    for (Integer msLevel : cache.keySet()){
      availableLevels.put(msLevel, true);
    }
    return availableLevels;
  }
  
  
  /**
   * calculates the shotgun intensity out of one chromatogram
   * @param cgChrom
   * @param shotgunType the type of shotgun processing (i.e. SHOTGUN_TYPE_MEAN, SHOTGUN_TYPE_MEDIAN, or SHOTGUN_TYPE_SUM)
   * @return calculated value
   */
  private static float calculateShotgunIntensity(CgChromatogram cgChrom, int shotgunType){
    //build the values of interest
    Vector<Float> ints = new Vector<Float>();
    float sum = 0f;
    for (int i=0; i!=cgChrom.Value.length; i++){
      sum += cgChrom.Value[i][1];
      if (!LipidomicsConstants.isShotgunIntensityRemoval() || cgChrom.Value[i][1]>0f)
        ints.add(cgChrom.Value[i][1]);
    }
    //if there is another threshold set, calculate the average, and remove the hits according to the cutoff
    if (LipidomicsConstants.isShotgunIntensityRemoval() && LipidomicsConstants.getShotgunRelIntCutoff()>0) {
      float cutoff = (sum*LipidomicsConstants.getShotgunRelIntCutoff())/((float)ints.size());
      Vector<Float> intermediate = new Vector<Float>();
      for (Float inten : ints) {
        if (inten>cutoff)
          intermediate.add(inten);
      }
    }
        
    float intensity = 0f;
    float[] intValues = new float[ints.size()];
    for (int i=0; i!=ints.size(); i++) intValues[i] = ints.get(i);
    if (shotgunType==SHOTGUN_TYPE_MEAN)
      intensity = Calculator.mean(intValues);
    else if (shotgunType==SHOTGUN_TYPE_MEDIAN)
      intensity = Calculator.median(intValues);
    else if (shotgunType==SHOTGUN_TYPE_SUM)
      intensity = sum;
    return intensity;
  }
  
}
