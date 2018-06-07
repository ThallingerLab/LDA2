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

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.ActionListener;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import at.tugraz.genome.dbutilities.Base64;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.LipidomicsChromatogram;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.lda.swing.RangeColor;
import at.tugraz.genome.lda.vos.SpectrumPointVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class Lipidomics2DSpectraChromPainter extends Lipidomics2DPainter
{
  private static final long serialVersionUID = 4970922671009881695L;
  
  private int pixelsForPeakDescription_ = 5;
  
  private Hashtable<Integer,String> spectra_;
  /** the MS2 precursor masses; key is the scan number; value is a vector of precorsor masses (higher MS levels have more than one precursor mass)*/
  private Hashtable<Integer,Vector<Double>> precursors_;
  /** lookup showing the MS-level for each scan number*/
  private Hashtable<Integer,Integer> levels_;
  private List<Integer> spectraSequence_;
  /** a list of scan numbers for MS2 spectra (higher MS-levels are not in this list)*/
  private List<Integer> ms2SpectraSequence_;
  
  private int spectrumSelected_;
  
  /** the default MS-level that is shown in the 3D-view*/
  private final static int DEFAULT_MS_LEVEL = 2;
  
  /** ranges where the spectrum has to be painted in a different color */
  private Hashtable<Integer,Vector<RangeColor>> rangeColors_;
  
  /** MS1 time ranges - checks if the colors shall be applied - if outside the range -> no coloring*/
  private Vector<Range> timeRanges_;
  
  private double annotationCutoff_;
  
  /** is this data shotgun data*/
  private boolean isShotgun_;
  
    
  /**
   * 
   * @param analyzer object for operations on the chrom data
   * @param rtNrSpectrumHash lookup table containing spectra in Base64 format; key is the scan number
   * @param rtNrPrecursorHash lookup table containing a vector of precursor masses for a scan; key is the scan number
   * @param scanNrLevelHash lookup table containing the MS-level the scan originates of; key is the scan number
   * @param retentionTimes lookup table containing the retention time of the scan; key is the scan number
   * @param precRt retention time of the MS-1 peak summit
   * @param mzStart for extracting chromatograms from raw data - here obsolete
   * @param mzStop for extracting chromatograms from raw data - here obsolete
   * @param resolutionFactor for extracting chromatograms from raw data - here obsolete
   * @param stepSize for extracting chromatograms from raw data - here obsolete
   * @param listener the listener who executes the commands
   * @param start for extracting chromatograms from raw data - here obsolete
   * @param stop for extracting chromatograms from raw data - here obsolete
   * @param raw applies for chromatograms - here obsolete
   * @param storedProbes applies for chromatograms - here obsolete
   * @param selectedProbes applies for chromatograms - here obsolete
   * @param isotopeNumber applies for chromatograms - here obsolete
   * @param charge applies for chromatograms - here obsolete
   * @param relativeIntensity show intensity in relative or absolute values
   * @param annotationCutoff the relative value where spectra peaks should not be annotated anymore
   * @param isShotgun is this shotgun data
   */
  public Lipidomics2DSpectraChromPainter(LipidomicsAnalyzer analyzer, Hashtable<Integer,String> rtNrSpectrumHash, Hashtable<Integer,Vector<Double>> rtNrPrecursorHash, 
      Hashtable<Integer, Integer> scanNrLevelHash, Hashtable<Integer,Float> retentionTimes, float precRt, float mzStart, float mzStop, int resolutionFactor,float stepSize, ActionListener listener,
      float start, float stop, boolean raw, Vector<CgProbe> storedProbes, Vector<CgProbe> selectedProbes,int isotopeNumber, int charge,
      boolean relativeIntensity, double annotationCutoff, boolean isShotgun){
    //TODO: The 2DSpectraChromPainter is currently not dedicated for storing areas - thus there is just 1 used as msLevel of the super constructor
    super(analyzer,new String[0], retentionTimes, mzStart, mzStop, resolutionFactor, stepSize, listener, DISPLAY_TIME_MZ,
        start, stop, raw, storedProbes, selectedProbes, isotopeNumber, charge,1,false);
    remove(m_txtMz);
    m_txtMz = null;
    spectra_ = rtNrSpectrumHash;
    precursors_ = rtNrPrecursorHash;
    levels_ = scanNrLevelHash;
    spectraSequence_ = new ArrayList<Integer>(spectra_.keySet());
    ms2SpectraSequence_ = new ArrayList<Integer>();
    Collections.sort(spectraSequence_);
    //only MS2 spectra are used for the sum spectrum
    for (Integer scanNr : spectraSequence_){
      if (scanNrLevelHash.get(scanNr)==2) ms2SpectraSequence_.add(scanNr);
    }
    spectrumSelected_ = -1;
    float diff = Float.MAX_VALUE;
    boolean containsDefaultLevel = false;
    for (Integer level : levels_.values()){
      if (level==DEFAULT_MS_LEVEL){
        containsDefaultLevel = true;
        break;
      }
    }
    for (int i=0; i!=spectraSequence_.size(); i++){
      int rtNr = spectraSequence_.get(i);
      if (containsDefaultLevel && levels_.get(rtNr)!=DEFAULT_MS_LEVEL)
        continue;
      float rt = retentionTimes_.get(rtNr);      
      if (Math.abs(rt-precRt)<diff){
        diff = Math.abs(rt-precRt);
        spectrumSelected_ = i;
      }
    }
    selectSpectrum(spectrumSelected_);
    m_minDispTime2d_ = mzStart;
    m_maxDispTime2d_ = mzStop;
    m_minTime_ = m_minDispTime2d_;
    m_maxTime_ = m_maxDispTime2d_;
    timeCorrectionFactor_ = 1f;
    relativeIntensity_ = relativeIntensity;
    rangeColors_ = null;
    timeRanges_ = null;
    m_2dGain_ = 0.97f;
    annotationCutoff_ = annotationCutoff;
    isShotgun_ = isShotgun;
//    spectrum_.LowerMzBand = this.stepSize_/2;
//    spectrum_.UpperMzBand = this.stepSize_/2;
//    spectrum_.Mz = (start+stop)/2;
//    spectrum_.Smooth(LipidomicsConstants.getChromSmoothRange(),
//        LipidomicsConstants.getChromSmoothRepeats());
//    spectrum_.GetMaximumAndAverage();
//
    
  }
  
  /**
   * constructor to create a Lipidomics2DSpectraChromPainter object - includes time ranges where spectrum is painted in a different color
   * @param analyzer analyzer to be used for peak calculation
   * @param rtNrSpectrumHash hash containing the spectra - key: scan number in chrom2 format; value: spectrum in Base64
   * @param rtNrPrecursorHash lookup table containing a vector of precursor masses for a scan; key is the scan number
   * @param scanNrLevelHash lookup table containing the MS-level the scan originates of; key is the scan number
   * @param retentionTimes retention time lookup - key is the scan number; value: retention time
   * @param precRt retention time of the precursor
   * @param mzStart for extracting chromatograms from raw data - here obsolete 
   * @param mzStop for extracting chromatograms from raw data - here obsolete
   * @param resolutionFactor for extracting chromatograms from raw data - here obsolete
   * @param stepSize for extracting chromatograms from raw data - here obsolete
   * @param listener the listener who executes the commands
   * @param start for extracting chromatograms from raw data - here obsolete
   * @param stop for extracting chromatograms from raw data - here obsolete
   * @param raw applies for chromatograms - here obsolete
   * @param storedProbes applies for chromatograms - here obsolete
   * @param selectedProbes applies for chromatograms - here obsolete
   * @param isotopeNumber applies for chromatograms - here obsolete
   * @param charge applies for chromatograms - here obsolete
   * @param relativeIntensity show intensity in relative or absolute values
   * @param param the detected LipidParameterSet to be displayed
   * @param rangeColors the m/z ranges that shall be painted in a different color
   * @param annotationCutoff the relative value where spectra peaks should not be annotated anymore
   * @param isShotgun is this shotgun data
   */
  public Lipidomics2DSpectraChromPainter(LipidomicsAnalyzer analyzer, Hashtable<Integer,String> rtNrSpectrumHash, Hashtable<Integer,Vector<Double>> rtNrPrecursorHash,
      Hashtable<Integer, Integer> scanNrLevelHash, Hashtable<Integer,Float> retentionTimes, float precRt, float mzStart, float mzStop, int resolutionFactor,float stepSize, ActionListener listener,
      float start, float stop, boolean raw, Vector<CgProbe> storedProbes, Vector<CgProbe> selectedProbes,int isotopeNumber, int charge,
      boolean relativeIntensity, LipidParameterSet param, Hashtable<Integer,Vector<RangeColor>> rangeColors, double annotationCutoff, boolean isShotgun){
    this(analyzer, rtNrSpectrumHash, rtNrPrecursorHash, scanNrLevelHash, retentionTimes, precRt, mzStart, mzStop, resolutionFactor, stepSize, listener,
        start, stop, raw, storedProbes, selectedProbes, isotopeNumber, charge, relativeIntensity, annotationCutoff, isShotgun);
    setRangeColors(param,rangeColors);
  }
  
  private void setRangeColors(LipidParameterSet param, Hashtable<Integer,Vector<RangeColor>> rangeColors){
    if (param instanceof LipidomicsMSnSet){
      rangeColors_ = rangeColors;
      timeRanges_ = new Vector<Range>();
      if (this.isShotgun_){
        timeRanges_.add(new Range(0f,Float.MAX_VALUE));
      }else{
        for (CgProbe probe : param.getIsotopicProbes().get(0)){
          timeRanges_.add(new Range(probe.LowerValley,probe.UpperValley));
        } 
      }
    }    
  }
  
  public void setParam(LipidParameterSet param)
  {
    if (param instanceof LipidomicsMSnSet)
    {      
      timeRanges_ = new Vector<Range>();
      if (this.isShotgun_){
        timeRanges_.add(new Range(0f,Float.MAX_VALUE));
      } else {
        for (CgProbe probe : param.getIsotopicProbes().get(0))
        {
          timeRanges_.add(new Range(probe.LowerValley,probe.UpperValley));
        }
      }
    }
  }
  
      
  protected void setInputLocations(){
    m_txtTime.setBounds(this.getWidth() - 100, 5, 95, 20);
    m_txtIntensity.setBounds(this.getWidth() - 100, 27, 95, 20);

    m_plusGain.setBounds(this.getWidth() - 65, m_txtIntensity.getY() + 23, 60,
        24);
    m_minusGain.setBounds(this.getWidth() - 65, m_plusGain.getY() + 25, 60, 24);
    m_plusShift.setBounds(this.getWidth() - 35, m_minusGain.getY() + 25, 30, 24);
    m_minusShift.setBounds(this.getWidth() - 65, m_minusGain.getY() + 25, 30,24);

  }
  
  @SuppressWarnings("unchecked")
  private LipidomicsChromatogram mergeSpectra(Collection<Integer> scanNumbers){
    Hashtable<String,SpectrumPointVO> intensityPoints = new Hashtable<String,SpectrumPointVO>();  
    for (Integer  scanNr : scanNumbers){
      String spectrum = this.spectra_.get(scanNr);
      FloatBuffer buffer = ByteBuffer.wrap(Base64.decode(spectrum)).asFloatBuffer();
      int iItems = buffer.limit();
      float mzValue = 0;
      for(int iItem = 0; iItem < iItems; iItem++){
        if (iItem%2==0) mzValue = buffer.get();
        else{
          float intensity = buffer.get();
          //TODO: I am not sure if I should use this reglementation (the if)
//            if (intensity>0){
          String mzOriginal = String.valueOf(mzValue);
          SpectrumPointVO vo = null;
          if (intensityPoints.containsKey(mzOriginal))
            vo = intensityPoints.get(mzOriginal);
          else
            vo = new SpectrumPointVO(mzOriginal,mzValue);
          vo.addIntensity(intensity);
          if (vo.getIntensity()> this.maxIntensity_)
            maxIntensity_ = vo.getIntensity();
          intensityPoints.put(mzOriginal, vo);
        }
//          }
      }
    }  
    List<SpectrumPointVO> vos = new ArrayList<SpectrumPointVO>(intensityPoints.values());
    Collections.sort(vos,new GeneralComparator("at.tugraz.genome.lda.vos.SpectrumPointVO", "getMz", "java.lang.Float"));
    LipidomicsChromatogram spectrum = new LipidomicsChromatogram(vos.size());
    for (int i=0;i!=vos.size();i++){
      SpectrumPointVO vo = vos.get(i);
      spectrum.Value[i][0] = vo.getMz();
      spectrum.Value[i][1] = vo.getIntensity();
    }
    return spectrum;
  }
  
  protected void draw2DDiagram(Graphics gx, float start, float stop, boolean raw){
    gx.setColor(Color.WHITE);
    gx.fillRect(0, 0, this.getWidth(), this.getHeight());
    gx.setColor(Color.BLACK);
//    chromMzStart_ = start;
//    chromMzStop_ = stop;
    raw_ = raw;
    int x0 = leftMargin_2d();
    int y0 = this.getHeight() - bottomMargin_2d();
    int w0 = this.diagramWidth_2d();
    int h0 = this.diagramHeight_2d();

//    Vector<CgProbe> paintableProbes = new Vector<CgProbe>();
//    for (int i = 0; i != this.storedProbes_.size(); i++) {
//      CgProbe aProbe = (CgProbe) this.storedProbes_.get(i);
//      boolean removed = false;
//      for (int j = 0; j != this.removedStoredProbes_.size(); j++) {
//        if (this.sameCgProbe((CgProbe) this.removedStoredProbes_.get(j), aProbe)) {
//          removed = true;
//        }
//      }
//      if (!removed)
//        paintableProbes.add(aProbe);
//    }
//
//
//    paint2dAreas(paintableProbes, Color.RED);
//    paint2dAreas(this.selectedProbes_, Color.GREEN);
    Vector<RangeColor> rColors = null;
    if (isSpectrumInIdentfiedRegion()){
      int msLevel = DEFAULT_MS_LEVEL;
      if (spectrumSelected_>-1)
        msLevel = levels_.get(spectraSequence_.get(spectrumSelected_));
      rColors = rangeColors_.get(msLevel);
    }
    draw2DDiagram(gx, cr_, x0,y0,w0,h0, m_minDispTime2d_, m_maxDispTime2d_,maxIntensity_,m_2dGain_,
        DISPLAY_TIME_MZ, 1, raw_,true,relativeIntensity_,true,null,rColors,false);
    drawSpectrumLegend(gx,cr_,0.5f,x0,y0,w0,h0,m_minDispTime2d_, m_maxDispTime2d_,maxIntensity_,m_2dGain_,
        raw,relativeIntensity_,rColors);
//    currentStart_ = startAsInt;
//    currentStop_ = stopAsInt;
//    this.m_txtMz.setText("mz="+Calculator.roundFloat((stop+start)/2f, 3)+"D");

  }
  
  /**
   * 
   * @return true if the spectrum is in a retention time range where a peak was found - m/z values can be colored 
   */
  private boolean isSpectrumInIdentfiedRegion(){
    if (spectrumSelected_<0) return true;
    else{
      if (this.timeRanges_!=null && timeRanges_.size()>0){
        float rt = retentionTimes_.get(spectraSequence_.get(spectrumSelected_));
        for (Range range : timeRanges_){
          if (range.insideRange(rt)) return true;
        }
      }
    }
    return false;
  }

  
  protected void drawSpectrumLegend(Graphics gx, LipidomicsChromatogram cr, float mzTolerance, 
      int x0, int y0, int w0, int h0, float mzStart, float mzStop, float maxIntensity, 
      float m_2dGain, boolean raw, boolean relativeValue,Vector<RangeColor> rColors){
    Vector<Integer> groupedInts = new Vector<Integer>();
    LinkedHashMap<Integer,Float> possIntsMz = new LinkedHashMap<Integer,Float>();
    LinkedHashMap<Integer,Float> possIntsInt = new LinkedHashMap<Integer,Float>();
    
    Hashtable<String,Float> highestInts = new Hashtable<String,Float>();
    Hashtable<String,Integer> highestIntPos = new Hashtable<String,Integer>();
    for (int i = 0; i < cr.Value.length; i++) {
      float mz = cr.Value[i][0]; 
      if (mz < mzStart)
        continue;
      if (mz > mzStop)
        continue;
      if (rColors!=null){
        RangeColor color = getColorAccordingToXPos(mz,rColors);
        if (color!=null){
          float highInt = 0f;
          if (highestInts.containsKey(color.getName())) highInt = highestInts.get(color.getName());
          float value = 0f;
          if (raw) value = cr.Value[i][1];
          else value = cr.Value[i][2];
          if (value>highInt){
            highInt = value;
            highestInts.put(color.getName(), highInt);
            highestIntPos.put(color.getName(), i);
          }
        }
      }
      for (Integer pos : possIntsMz.keySet()){
        // The value is outside the tolerance value and does not need to be grouped
        if ((possIntsMz.get(pos)+mzTolerance)<mz){
          groupedInts.add(pos);
          possIntsMz.remove(pos);
          possIntsInt.remove(pos);
        }
      }
      float value = 0f;
      if (raw) value = cr.Value[i][1];
      else value = cr.Value[i][2];
      if (value>0){
        boolean foundHigherOne = false;
        for (Integer pos : possIntsMz.keySet()){
          float otherValue = possIntsInt.get(pos);
          if (otherValue>value)
            foundHigherOne = true;
          else{
            possIntsMz.remove(pos);
            possIntsInt.remove(pos);
          }
        }

        if (!foundHigherOne){
          possIntsMz.put(i,mz);
          possIntsInt.put(i,value);
        }
      }
    }
    for (Integer pos : possIntsMz.keySet())groupedInts.add(pos);
    gx.setFont(new Font("SansSerif", Font.PLAIN, 10));
    FontMetrics fm = gx.getFontMetrics();
    for (Integer pos : groupedInts){
      int[] coords = getCoordinatesInDiagram(cr.Value, pos, x0, y0, w0, h0,
          mzStart, mzStop, maxIntensity, m_2dGain, raw, relativeValue,false);
      // these ones should get a description
      double rel = ((double)(cr.Value[pos][1]/maxIntensity))*100d;
      //check if it is a found hit
      boolean isAnnotatedPeak = false;
      RangeColor color = null;
      if (rColors!=null){
        for (RangeColor rc : rColors){
          if (!highestIntPos.containsKey(rc.getName())) continue;
           float rt = cr.Value[pos][0];
          float highestRt = cr.Value[highestIntPos.get(rc.getName())][0];
          if ((rt-mzTolerance)<highestRt && highestRt<(rt+mzTolerance)){
            isAnnotatedPeak=true; 
            color = rc;
          }
        }
      }
      if (coords[0]>x0&&coords[0]<(x0+w0) && coords[1]>(y0-h0)&&coords[1]<(y0-h0/500) && (rel>=annotationCutoff_ || isAnnotatedPeak)){
        String mzValue = String.valueOf(Calculator.roundFloat(cr.Value[pos][0], 3));
        int x = coords[0]-fm.stringWidth(mzValue)/2;
        int y = coords[1]- pixelsForPeakDescription_-fm.getDescent();
        if (isAnnotatedPeak) gx.setColor(color.getColor());
        gx.drawString(mzValue, x, y);
        if (isAnnotatedPeak) gx.setColor(Color.BLACK);
      }  
    }
    if (rColors!=null){
      for (RangeColor rc : rColors){
        if (!highestIntPos.containsKey(rc.getName())) continue;
        int[] coords = getCoordinatesInDiagram(cr.Value, highestIntPos.get(rc.getName()), x0, y0, w0, h0,
            mzStart, mzStop, maxIntensity, m_2dGain, raw, relativeValue,false);
        // these ones should get a description
        if (coords[0]>x0&&coords[0]<(x0+w0)&&coords[1]>(y0-h0)&&coords[1]<(y0-h0/500)){
          int x = coords[0]-fm.stringWidth(rc.getName())/2;
          int y = coords[1]-fm.getHeight()-pixelsForPeakDescription_-fm.getDescent();

          gx.setColor(rc.getColor());
          gx.drawString(rc.getName(), x, y);
        }  
      }
      gx.setColor(Color.BLACK);
    }
  }

  public void setShownRetentionTime(float newVal)
  {
    m_txtTime.setText("m/z = " + Calculator.roundFloat(newVal,4));
  }
  
  public int checkPopupClick(int x, int y){
    return 0;
  }
  
  public void nextSpectrum(LipidParameterSet param, Hashtable<Integer,Vector<RangeColor>> rangeColors){
    setRangeColors(param,rangeColors);
    selectSpectrum(this.spectrumSelected_+1);
  }

  public void nextSpectrum(){
    selectSpectrum(this.spectrumSelected_+1);
  }
  
  public void previousSpectrum(LipidParameterSet param, Hashtable<Integer,Vector<RangeColor>> rangeColors){
    setRangeColors(param,rangeColors);
    previousSpectrum();
  }

  public void refresh(LipidParameterSet param, Hashtable<Integer,Vector<RangeColor>> rangeColors){
    setRangeColors(param,rangeColors);
    selectSpectrum(this.spectrumSelected_);
  }
  
  /**
   * causes the colored peaks to be cleared, even if
   * the new LipidParameterSet is not of type LipidomicsMSnSet
   */
  public void clearRangeColors(){
    rangeColors_ = new Hashtable<Integer,Vector<RangeColor>>();
    timeRanges_ = new Vector<Range>();
    selectSpectrum(this.spectrumSelected_);
  }
  
  public void previousSpectrum(){
    selectSpectrum(this.spectrumSelected_-1);    
  }
  
  private void selectSpectrum(int spectrumNr){
    spectrumSelected_ = spectrumNr;
    if (spectrumSelected_>=this.spectra_.size()) this.spectrumSelected_ = -1;
    if (spectrumSelected_<-1) this.spectrumSelected_ = this.spectra_.size()-1;
    maxIntensity_ = 0;
    if (spectrumSelected_==-1){
      cr_ = mergeSpectra(this.ms2SpectraSequence_);
    } else{
      Vector<Integer> scanNumbers = new Vector<Integer>();
      int scanNr = spectraSequence_.get(spectrumSelected_);
      scanNumbers.add(scanNr);
      cr_ = mergeSpectra(scanNumbers);
    }
    this.repaint();
  }
  
  public String getSpectSelectedText(){
    if (spectrumSelected_==-1) return "Sum("+String.valueOf(this.ms2SpectraSequence_.size())+")";
    else return String.valueOf(this.spectrumSelected_+1)+"/"+this.spectra_.size();
  }
  
  public String getRtSelectedText(){
    if (spectrumSelected_==-1) return "";
    else{
      return String.valueOf(Calculator.roundFloat(retentionTimes_.get(spectraSequence_.get(spectrumSelected_))/60f,2));
    }
  }
  
  public Vector<Double> getPrecursorMassSelected(){
    if (spectrumSelected_==-1) return new Vector<Double>();
    else{
      return precursors_.get(spectraSequence_.get(spectrumSelected_));
    }
  }

  public String getMsLevelSelected(){
    if (spectrumSelected_==-1) return "";
    else{
      return String.valueOf(levels_.get(spectraSequence_.get(spectrumSelected_)));
    }
  }

  
  public float[] getRTRange(){
    if (spectrumSelected_==-1)return(new float[]{0f,0f});
    float[] range = new float[2];
    int startSpectrumNr = spectrumSelected_;
    while (startSpectrumNr!=0 && levels_.get(spectraSequence_.get(startSpectrumNr))!=DEFAULT_MS_LEVEL)
      startSpectrumNr--;
    if (startSpectrumNr<1) range[0] = retentionTimes_.get(spectraSequence_.get(0));
    else{
      startSpectrumNr--;
      while (startSpectrumNr!=0 && levels_.get(spectraSequence_.get(startSpectrumNr))!=DEFAULT_MS_LEVEL)
        startSpectrumNr--;
      range[0] = retentionTimes_.get(spectraSequence_.get(startSpectrumNr));
    }
    int stopSpectrumNr = spectrumSelected_;
    while (stopSpectrumNr!=(spectraSequence_.size()-1) && levels_.get(spectraSequence_.get(stopSpectrumNr))!=DEFAULT_MS_LEVEL)
      stopSpectrumNr++;
    if (stopSpectrumNr>=(spectraSequence_.size()-1)) range[1] = retentionTimes_.get(spectraSequence_.get(spectraSequence_.size()-1));
    else{
      stopSpectrumNr++;
      while (stopSpectrumNr!=(spectraSequence_.size()-1) && levels_.get(spectraSequence_.get(stopSpectrumNr))!=DEFAULT_MS_LEVEL)
        stopSpectrumNr++;
      range[1] = retentionTimes_.get(spectraSequence_.get(stopSpectrumNr));
    }    
    return range;
  }
  
  public int getSpectrumSelected(){
    return this.spectrumSelected_;
  }
  
  public int getMs2LevelSpectrumSelected(){
    return getMs2LevelSpectrumSelected(spectrumSelected_);
  }
  
  public int getMs2LevelSpectrumSelected(int specNr){
    int ms2Number = -1;
    int nr = specNr;
    if (nr>=this.spectra_.size()) nr = -1;
    if (nr<-1) nr = this.spectra_.size()-1;
    for (int i=0; i!=(nr+1); i++){
      if (levels_.get(spectraSequence_.get(i))!=DEFAULT_MS_LEVEL)
        continue;
      ms2Number++;
    }
    return ms2Number;
  }
  
  public void setAnnotationThreshold(double cutoff){
    this.annotationCutoff_ = cutoff;
    selectSpectrum(spectrumSelected_);
  }


/*  protected static void drawDiagramItself(Graphics gx, CgChromatogram cr, int x0, int y0, int w0, int h0, float m_minDispTime2d, 
      float m_maxDispTime2d, float maxIntensity, float m_2dGain, boolean raw, boolean relativeValue){
    int x = 0, y = 0;
    gx.setColor(Color.BLACK);
    boolean firstPoint = true;
    int len = cr.Value.length;
    for (int i = 0; i < len; i++) {

      if (cr.Value[i][0] < m_minDispTime2d)
        continue;
      if (cr.Value[i][0] > m_maxDispTime2d)
        continue;
      int[] coords = getCoordinatesInDiagram(cr.Value, i, x0, y0, w0, h0, m_minDispTime2d, 
          m_maxDispTime2d, maxIntensity, m_2dGain, raw, relativeValue);
      
      x = coords[0];
      y = coords[1];
      if (firstPoint == true) {
        firstPoint = false;
      }
      gx.drawLine(x, h0, x, y);
    }

  }*/

  public void writeMgf(File file, String rawFileName) throws IOException{
    BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(file.getAbsolutePath()));
    stream.write(("BEGIN IONS\n").getBytes());
    stream.write(("TITLE="+rawFileName+"\n").getBytes());
    String rtString = "RTINSECONDS=";
    String pepMassString = "PEPMASS=";
    if (spectrumSelected_>=0){
      rtString += Calculator.roundFloat(retentionTimes_.get(spectraSequence_.get(spectrumSelected_)),4)+"\n";
      pepMassString += precursors_.get(spectraSequence_.get(spectrumSelected_))+"\n";
    }else{
      rtString += "\n";
      pepMassString += "\n";
      
    }
    stream.write(rtString.getBytes());
    stream.write(pepMassString.getBytes());
    for (int i=0; i!=cr_.Value.length; i++){
      if (cr_.Value[i][1]>0) stream.write((cr_.Value[i][0]+" "+cr_.Value[i][1]+"\n").getBytes());
    }
    stream.write(("END IONS\n").getBytes());
    stream.close();
  }
  
}
