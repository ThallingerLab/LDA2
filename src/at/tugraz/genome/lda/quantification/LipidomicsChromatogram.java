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

package at.tugraz.genome.lda.quantification;

import java.util.Collections;
import java.util.Vector;

import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.util.FloatMatrix;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsChromatogram extends CgChromatogram
{
  
  private boolean greedyFragment_;
  private int loSteepnessPoint_;
  private int upSteepnessPoint_;
  
  // the steepness changes for the peak borders
  public float greedySteepnessChange1_;
  public float greedySteepnessChange2_;
  
  // there are 2 relative intensity values where different steepness
  // changes are regarded as peak border
  public float greedyIntensityCutoff1_;
  public float greedyIntensityCutoff2_;

  // for the peak detection a cut-off value relative to the highest
  // value can be specified to set the peak border
  // just valid for GreedySteepnessImproved method
  public float intensityCutoff_;

  
  public LipidomicsChromatogram(int size){
    super(size);
    greedyFragment_=false;
    greedySteepnessChange1_ = 1.5f;
    greedyIntensityCutoff1_ = 0.15f;
    greedySteepnessChange2_ = 1.8f;
    greedyIntensityCutoff2_ = 0.2f;
    intensityCutoff_ = 0.00f;
  }
  
  public LipidomicsChromatogram(CgChromatogram cgc){
    this(cgc.Value.length);
    LowerMzBand = cgc.LowerMzBand;// Lower m/z area
    Mz = cgc.Mz;// m/z value
    UpperMzBand = cgc.UpperMzBand;// Upper m/z area
    ScanCount = cgc.ScanCount;// Number of scans
    Value = cgc.Value.clone();// The Chromatogram...
    LoValley = cgc.LoValley;// Index of lower valley
    Peak = cgc.Peak;// Index of peak;
    UpValley = cgc.UpValley;// Index of upper valley
    PeakTime = cgc.PeakTime;// Retention Time of Peak
    PeakTimeAverage = cgc.PeakTimeAverage;// Retention Time of Peak averaged
    Background = cgc.Background;// The background
    Area = cgc.Area;// The area of our peak
    AreaErr = cgc.AreaErr;// The error of the area
    Good = cgc.Good;// Does the area make sense
    m_peakValue = cgc.getM_peakValue();// Chroma peak (internal)
    m_average = cgc.getM_average();// Chroma average (internal)
    startPosition_ = cgc.getStartPosition();
    isProfile_ = cgc.isProfile_;
    highestIntensity_ = cgc.getHighestIntensity();
    startSmoothScan_ = cgc.startSmoothScan_;
    stopSmoothScan_ = cgc.stopSmoothScan_;
  }
  
  public void FindPeak(int position, int direction, int valleyMethod)
  {
    this.greedyFragment_ = false;
    loSteepnessPoint_ = -1;
    upSteepnessPoint_ = -1;
    super.FindPeak(position, direction, valleyMethod);
  }
  
  protected void GetValleys(int method)
  {
    if (method==LipidomicsDefines.StandardValleyMethod)
          GetValleysOriginal();
      else if (method==LipidomicsDefines.EnhancedValleyMethod) 
          GetValleysEnhanced(true);
      else if (method==LipidomicsDefines.GreedySteepnessReductionMethod) 
        GetValleysGreedy();
      else if (method==LipidomicsDefines.GreedySteepnessPoint) 
        GetValleysGreedySteepnessPoint();
      else if (method==LipidomicsDefines.GreedySteepnessImproved) 
        GetValleysGreedyImproved();
      else if (method==LipidomicsDefines.GreedySteepnessPointImproved) 
        GetValleysGreedySteepnessPointImproved();
      else
          GetValleysOriginal();
  }
  
  public void Smooth(float range, int repeats, boolean copyRawDataFirst, SavGolJNI sav_gol)
  {
    if (copyRawDataFirst) copyRawData();

    float result[];
    result = sav_gol.Smooth(  Value, ScanCount, range, repeats, startSmoothScan_, stopSmoothScan_ );
    for( int i = 0; i < ScanCount; i ++ )
    {
      Value[i][2] = result[i];
    }
  }

  public void Smooth(float range, int repeats, SavGolJNI sav_gol)
  {
    this.Smooth(range, repeats, true, sav_gol);
  }
  
  /**
   * This function finds the next valley from left (when "direction 
   * = -1") or right (when "direction = 1"), starting at "position" 
   * in the smoothed chromatogram. It accounts for quick changes in 
   * the steepness of a peek; then a potential overlap must be assumed
   */
  private void GetValleysGreedy()
  {
    GetValleysGreedySteepnessPoint();
    int steepnessPointLoValley = LoValley;
    int steepnessPointUpValley = UpValley;
      // if I detect an overlap then I have to set the border of the one peak exactly in the midst to the other peak
      // I do this by detection the potential peak maximum of the overlapping peak, by detecting when the peak starts
      // to fall again steeper
    int position = UpValley;
    float[] diffs = this.calculateUpwardsDifferential(position);
    float differential1 = diffs[0];
    float differential2 = diffs[1];
    if ((1.5*differential2)<differential1){
      int oldPosition = position;
      while(position<ScanCount-1)
      {
        if (position>0)
          differential1 = Value[position-1][2]-Value[position][2];
          differential2 = Value[position][2]-Value[position+1][2];
          if (differential2>(1.1*differential1)||Value[position+1][2]>Value[position][2])break;
          position++;
        }
      // if the peak intensity of the presumeable peak is bigger 10% then we can assume that 
      // we have found a different peak; then take the midst between peak an steepness reduction point;
      // otherwise the steepness reduction point is the valley
      if (Value[Peak][2]*0.1<Value[position][2]){
        position = (position+oldPosition)/2+1;
      }else
        position = oldPosition;
    }
    UpValley = position;
    if (UpValley>ScanCount-2) UpValley = ScanCount-2;

    // if I detect an overlap then I have to set the border of the one peak exactly in the midst to the other peak
    // I do this by detection the potential peak maximum of the overlapping peak, by detecting when the peak starts
    // to fall again steeper
    position = LoValley;
    diffs = this.calculateDownwardsDifferential(position);
    differential1 = diffs[0];
    differential2 = diffs[1];
    if ((1.5*differential2)<differential1){
      int oldPosition = position;
      while(position>0)
      {
        differential1 = Value[position+1][2]-Value[position][2];
        differential2 = Value[position][2]-Value[position-1][2];
        if (differential2>(1.1*differential1)||Value[position-1][2]>Value[position][2])break;
        position--;
      }
      // if the peak intensity of the presumeable peak is bigger 10% then we can assume that 
      // we have found a different peak; then take the midst between peak an steepness reduction point;
      // otherwise the steepness reduction point is the valley
        if (Value[Peak][2]*0.1<Value[position][2]){
          position = (position+oldPosition)/2;
        }else{
          position = oldPosition;
        }  
      }

    //check if the peak we are looking at is a small adjacent peak
      LoValley = position;
      if (LoValley<1) LoValley = 1;
      if (this.startPosition_<steepnessPointLoValley || steepnessPointUpValley<this.startPosition_){
        int greedyLoValley = LoValley;
        int greedyUpValley = UpValley;
        int peakDistanceGreedy = greedyUpValley-greedyLoValley;
        this.GetValleysOriginal();
        int distanceRest = 0;
        if (this.startPosition_<steepnessPointLoValley){
          distanceRest = greedyLoValley - LoValley;
        }  
        if (this.startPosition_>steepnessPointUpValley){
          distanceRest = UpValley - greedyUpValley;
        }
        // if the peak is big enough, put it in the rest between greedy peak and remaining part
        if (distanceRest>(peakDistanceGreedy)/3){
          this.greedyFragment_ = true;
          if (this.startPosition_<steepnessPointLoValley){
            UpValley = greedyLoValley;             
          }
          if (this.startPosition_>steepnessPointUpValley){
            LoValley = greedyUpValley;
          }
          this.Peak = (UpValley+LoValley)/2;
        }else{
          LoValley = greedyLoValley;
          UpValley = greedyUpValley;
        }
      }
  }
  
  private void GetValleysGreedySteepnessPoint(){
    int position;
    // =============================================================
    // Do the upwards direction:
    // =============================================================
    position = Peak;
    float differential1 = 0;
    float differential2 = 0;
    while(position<ScanCount-1)
    {
      float[] diffs = this.calculateUpwardsDifferential(position);
      differential1 = diffs[0];
      differential2 = diffs[1];
      if (position>=ScanCount-2){
        position++;
        break;
      }  
      if ((1.5*differential2)<differential1||Value[position+1][2]>Value[position][2]) break;
      position++;
    }
    UpValley = position;
    if (UpValley>ScanCount-2) UpValley = ScanCount-2;

    // =============================================================
    // Do the downwards direction:
    // =============================================================
    
    position = Peak;
    differential1 = 0;
    differential2 = 0;
    while(position>0)
    {
      
      float[] diffs = this.calculateDownwardsDifferential(position);
      differential1 = diffs[0];
      differential2 = diffs[1];
      if (position<ScanCount-2)
        differential1 = Value[position+1][2]-Value[position][2];
      if (position>0)
        differential2 = Value[position][2]-Value[position-1][2];
      else{
        break;
      }

      if ((1.5*differential2)<differential1||Value[position-1][2]>Value[position][2]) break;
      position--;
    }
    LoValley = position;
    if (LoValley<1) LoValley = 1;
  }
    
  private void GetValleysGreedyImproved(){
    this.GetValleysGreedySteepnessPointImproved();
    int steepnessPointLoValley = LoValley;
    int steepnessPointUpValley = UpValley;
    //check if the peak we are looking at is a small adjacent peak
    if (LoValley<1) LoValley = 1;
    if ((this.startPosition_<steepnessPointLoValley&&steepnessPointLoValley!=-1)|| (steepnessPointUpValley<this.startPosition_&&steepnessPointUpValley!=-1)){
      int greedyLoValley = LoValley;
      int greedyUpValley = UpValley;
      int peakDistanceGreedy = greedyUpValley-greedyLoValley;
      this.GetValleysOriginal();
      int distanceRest = 0;
      if (this.startPosition_<steepnessPointLoValley){
        distanceRest = greedyLoValley - LoValley;
      }  
      if (this.startPosition_>steepnessPointUpValley){
        distanceRest = UpValley - greedyUpValley;
      }
      // if the peak is big enough, put it in the rest between greedy peak and remaining part
      if (distanceRest>(peakDistanceGreedy)/3){
        this.greedyFragment_ = true;
        if (this.startPosition_<steepnessPointLoValley){
          UpValley = greedyLoValley;
          if (this.loSteepnessPoint_!=-1&&this.loSteepnessPoint_>UpValley)
            this.UpValley = loSteepnessPoint_+(loSteepnessPoint_-UpValley);
        }
        if (this.startPosition_>steepnessPointUpValley){
          LoValley = greedyUpValley;
          if (this.upSteepnessPoint_!=-1&&this.upSteepnessPoint_<LoValley)
            this.LoValley = upSteepnessPoint_+(LoValley-upSteepnessPoint_);
        }
        this.Peak = (UpValley+LoValley)/2;
        if (this.startPosition_<steepnessPointLoValley){
          this.upSteepnessPoint_ = this.loSteepnessPoint_;
          this.loSteepnessPoint_ = -1;
        } 
        if (this.startPosition_>steepnessPointUpValley){
          this.loSteepnessPoint_ = this.upSteepnessPoint_;
          this.upSteepnessPoint_ = -1;
        }
      }else{
        LoValley = greedyLoValley;
        UpValley = greedyUpValley;
      }
    }
    if (UpValley>ScanCount-2) UpValley = ScanCount-2;
    if (LoValley<1) LoValley = 1;
  }
  
  private void GetValleysGreedySteepnessPointImproved(){
    int position;
    // =============================================================
    // Do the upwards direction:
    // =============================================================
    position = Peak;
    // only a security check if nothing is there
    if (Value[position][2]<0.000000001f){
      LoValley = position-1;
      UpValley = position+1;
      if (UpValley>ScanCount-2){
        UpValley = ScanCount-2;
        LoValley = ScanCount-4;
      }
      if (LoValley<1){
        LoValley = 1;
        UpValley = 3;
      }
      return;
    }
    float differential1 = 0;
    float differential2 = 0;
//    float highestDifferential = 0f;
    // this method starts at the position of the peak and goes upwards
    // for every position a differential quotient is calculated
    // if the differential quotient changes too quickly (depending on the relative intensity)
    // an overlapping peak border is assumed
    // if the value is lower than a specific cutoff than a border is assumed
    // if the value of the next point is higher than of the previous one a border is assumed      
    while(position<ScanCount-1)
    {
      float[] diffs = this.calculateUpwardsDifferential(position);
      differential1 = diffs[0];
      differential2 = diffs[1];
//      if (highestDifferential == 0f)
//        differential1 = highestDifferential;
      if (position>=ScanCount-2){
        position++;
        break;
      }
      //ignore changes that are lower than 1 per mille - for peak plateaus
      if (Math.abs(differential1)<0.001f*Value[position][2]){
        position++;
        continue;
      }
      // for a border position:
      // below greedyIntensityCutoff1_ no steepness change is regarded just if the next value is higher than the previous one
      // between  greedyIntensityCutoff1_ and greedyIntensityCutoff2_ the steepness change has to be higher than greedySteepnessChange2_
      // the reason: in smaller regions higher steepness changes are possible without being a peak border
      // over greedyIntensityCutoff2_ the the steepness change has to be higher than greedySteepnessChange1_
      if (((greedySteepnessChange1_*differential2)<differential1&&(Value[position][2]/Value[Peak][2])>greedyIntensityCutoff1_)&&
          ((greedySteepnessChange2_*differential2)<differential1||(Value[position][2]/Value[Peak][2])>greedyIntensityCutoff2_)){
        // stores the point of the steepness change
        this.upSteepnessPoint_ = position;
        // calculates the a line prolongating the current steepness
        // where the line hits the 0 intensity line the peak border is assumed
        if (differential1>0){
          float differential = differential1;
//          if (position>2 && (Value[position-2][2]-Value[position-1][2])>differential)
//            differential =  Value[position-2][2]-Value[position-1][2];
//          if (position>3 && (Value[position-3][2]-Value[position-2][2])>differential)
//            differential =  Value[position-3][2]-Value[position-2][2];
          float newDifferential = getHighestCloseDifferentialUp(position);
          if (newDifferential>differential)
            differential = newDifferential;
          int positionDifference = (int)Calculator.roundFloat((Value[position][2])/differential,0);
          position = position+positionDifference;
        }
        // this is necessary if there is an overlap to an higher peak and the steepness is quite low
        // if this would not be done a lot of the higher peaks would be taken to the smaller one
        int interimPosition = this.upSteepnessPoint_+1;
        while (interimPosition<position&&interimPosition<(ScanCount-1)){
          if (Value[interimPosition+1][2]>Value[interimPosition][2]||Value[interimPosition][2]<1){
            position = interimPosition;
            break;
          }
          interimPosition++;
        }
        if (position>(ScanCount-1))
          position = (ScanCount-1);
        break;
      // if the next value is higher than the previous one or the value is below the intensity-cutoff a border is assumed
      }else if (Value[position+1][2]>Value[position][2]||Value[position+1][2]<Value[Peak][2]*intensityCutoff_/*||differential1<(highestDifferential/1000)*/)
        break;
//      if (differential1>highestDifferential)
//        highestDifferential = differential1;
      position++;
    }     
    UpValley = position;
    if (UpValley>ScanCount-2) UpValley = ScanCount-2;
    // =============================================================
    // Do the downwards direction:
    // =============================================================
    
    position = Peak;
    differential1 = 0;
    differential2 = 0;
//    highestDifferential = 0f;
    // this method starts at the position of the peak and goes upwards
    // for every position a differential quotient is calculated
    // if the differential quotient changes too quickly (depending on the relative intensity)
    // an overlapping peak border is assumed
    // if the value is lower than a specific cutoff than a border is assumed
    // if the value of the next point is higher than of the previous one a border is assumed
    while(position>0)
    {
      
      float[] diffs = this.calculateDownwardsDifferential(position);
      differential1 = diffs[0];
      differential2 = diffs[1];
//      if (highestDifferential == 0f)
//        differential1 = highestDifferential;
      if (position<ScanCount-2)
        differential1 = Value[position+1][2]-Value[position][2];
      if (position>0)
        differential2 = Value[position][2]-Value[position-1][2];
      else{
        break;
      }
      //ignore changes that are lower than 1 per mille - for peak plateaus
      if (Math.abs(differential1)<0.001f*Value[position][2]){
        position--;
        continue;
      }
      // for a border position:
      // below greedyIntensityCutoff1_ no steepness change is regarded just if the next value is higher than the previous one
      // between  greedyIntensityCutoff1_ and greedyIntensityCutoff2_ the steepness change has to be higher than greedySteepnessChange2_
      // the reason: in smaller regions higher steepness changes are possible without being a peak border
      // over greedyIntensityCutoff2_ the the steepness change has to be higher than greedySteepnessChange1_
      if (((greedySteepnessChange1_*differential2)<differential1&&(Value[position][2]/Value[Peak][2])>greedyIntensityCutoff1_)&&
          ((greedySteepnessChange2_*differential2)<differential1||(Value[position][2]/Value[Peak][2])>greedyIntensityCutoff2_)){
        // stores the point of the steepness change
        this.loSteepnessPoint_ = position;
        // calculates the a line prolongating the current steepness
        // where the line hits the 0 intensity line the peak border is assumed
        if (differential1>0){
          float differential = differential1;
          float newDifferential = getHighestCloseDifferentialLo(position);
          if (newDifferential>differential)
            differential = newDifferential;
          int positionDifference = (int)Calculator.roundFloat((Value[position][2])/differential,0);
          position = position-positionDifference;
        }
        // this is necessary if there is an overlap to an higher peak and the steepness is quite low
        // if this would not be done a lot of the higher peak would be taken to the smaller one
        int interimPosition = this.loSteepnessPoint_-1;
        while (interimPosition>position&&interimPosition>0){
          if (Value[interimPosition-1][2]>Value[interimPosition][2]||Value[interimPosition][2]<1){
            position = interimPosition;
            break;
          }
          interimPosition--;
        }
        if (position<0)
          position = 0;
        break;            
      }else if (Value[position-1][2]>Value[position][2]||Value[position-1][2]<Value[Peak][2]*intensityCutoff_/*||differential1<(highestDifferential/100)*/){
//        for (int i=0;i!=20;i++){
//          System.out.println(Value[position-i][0]+";"+Value[position-i][2]);
//        }
        break;
      }
//      if (differential1>highestDifferential)
//        highestDifferential = differential1;
      position--;
    }
    LoValley = position;
    if (LoValley<1) LoValley = 1;
    // this is just for prevention of errors
    if  (((Peak-LoValley)>((UpValley-Peak)*7)&&(LoValley<Peak)&&(Peak<UpValley))||((UpValley-Peak)>((Peak-LoValley)*7)&&(LoValley<Peak)&&(Peak<UpValley))){
      this.GetBackground(50);
      if ((Peak-LoValley)>((UpValley-Peak)*7)&&(LoValley<Peak)&&(Peak<UpValley)){
        LoValley = Peak-(UpValley-Peak);
        while (LoValley>-1 && Value[LoValley][2]>Value[Peak][2]*0.05 && Value[LoValley][2]>this.Background*3 ) LoValley--; 
      }
      if ((UpValley-Peak)>((Peak-LoValley)*7)&&(LoValley<Peak)&&(Peak<UpValley)){
        UpValley = Peak+(Peak-LoValley);
        while (UpValley<ScanCount && Value[UpValley][2]>Value[Peak][2]*0.05 && Value[UpValley][2]>this.Background*3 ) UpValley++;
      }
    }
    if (UpValley>ScanCount-2) UpValley = ScanCount-2;
    if (LoValley<1) LoValley = 1;
  }
  
  /**
   * detects the most intense signal inside a given chromatogram border 
   * @param low lower border
   * @param high upper border
   */
  public void detectPeakIntensityInsideBorders(int low, int high){
    LoValley = low;
    UpValley = high;
    float intensity = 0f;
    for (int i=LoValley; i<UpValley; i++){
      if (i<0||i>=ScanCount) continue;
      if (Value[i][2]>intensity){
        intensity = Value[i][2];
        Peak = i;
      }
    }
  }

  
  private float getHighestCloseDifferentialUp(int position){
    float differential = 0f;
    for (int i=0;i!=10;i++){
      if (position-i-1>0 && (Value[position-i-1][2]-Value[position-i][2])>differential)
        differential =  Value[position-i-1][2]-Value[position-i][2];
    }
    return differential;
  }
  
  
  private float getHighestCloseDifferentialLo(int position){
    float differential = 0f;
    for (int i=0;i!=10;i++){
      if (position<ScanCount-(i+1) && (Value[position+i+1][2]-Value[position+i][2])>differential)
        differential =  Value[position+i+1][2]-Value[position+i][2];
    }
    return differential;
  }
  
  private float[] calculateDownwardsDifferential(int position){
    float[] differentials = new float[2];
    float differential1 = 0;
    float differential2 = 0;
    if (position<ScanCount-2)
      differential1 = Value[position+1][2]-Value[position][2];
    if (position>0)
      differential2 = Value[position][2]-Value[position-1][2];
    differentials[0] = differential1;
    differentials[1] = differential2;
    return differentials;
  }
  
  private float[] calculateUpwardsDifferential(int position){
    float[] differentials = new float[2];
    float differential1 = 0;
    float differential2 = 0;
    if (position>0)
      differential1 = Value[position-1][2]-Value[position][2];
    if (position<ScanCount-2)
      differential2 = Value[position][2]-Value[position+1][2];
    differentials[0] = differential1;
    differentials[1] = differential2;
    return differentials;
  }
  
  public boolean isGreedyFragment()
  {
    return greedyFragment_;
  }


  public int getLoSteepnessPoint()
  {
    return loSteepnessPoint_;
  }


  public int getUpSteepnessPoint()
  {
    return upSteepnessPoint_;
  }

  /**
   * returns peak borders, where the intensities of the peak drop below a certain percentual threshold
   * @param percent the percentage for the threshold
   * @param start the start of the peak (conventional borders)
   * @param stop the end of the peak (conventional borders)
   * @return float[0] = border start value; float[1] = border stop value
   */
  public float[] findPercentualBorders(float percent, float start, float stop){
    float startPercent = Float.MAX_VALUE;
    float stopPercent = 0f;
    int startScan = LipidomicsAnalyzer.findIndexByTime(start, this);
    startScan--;
    if (startScan<0) startScan++; 
    int stopScan = LipidomicsAnalyzer.findIndexByTime(stop, this);
    stopScan++;
    if (stopScan>=this.Value.length) stopScan--;
    int highestScanNumber = -1;
    float intensity = 0f;
    for (int i=startScan; i<=stopScan; i++){
      if (Value[i][2]>intensity){
        highestScanNumber = i;
        intensity = Value[i][2];
      }
    }
    
    float threshold = (intensity*percent)/100f;
    if (highestScanNumber>-1){
      //scan in negative time direction
      int scanCount = highestScanNumber;
      for (scanCount=highestScanNumber; scanCount>=startScan; scanCount--){
        if (Value[scanCount][2]<threshold) break;
      }
      if (scanCount<=startScan) startPercent = start;
      else if (scanCount<0) startPercent = 0f;
      else if ((scanCount+1)==Value.length) startPercent = Value[scanCount][0];
      else {
        startPercent = Value[scanCount][0] + (Value[scanCount+1][0]-Value[scanCount][0])/10f; 
      }
      //scan in positive time direction
      scanCount = highestScanNumber;
      for (scanCount=highestScanNumber; scanCount<=stopScan; scanCount++){
        if (Value[scanCount][2]<threshold) break;
      }
      if (scanCount>=stopScan) stopPercent = stop;
      else if (scanCount>=Value.length-1) stopPercent = Float.MAX_VALUE;
      else if (scanCount<=0) stopPercent = Value[0][0];
      else {
        stopPercent = Value[scanCount][0] +(Value[scanCount][0]-Value[scanCount-1][0])/10f;
      }
    }
    float[] result = new float[2];
    result[0] = startPercent;
    result[1] = stopPercent;
    return result;
  }
  
  
  /**
   * 
   * @return the index of the highest intensity in a chromatogram
   */
  public int findIndexOfHighestIntensity(){
    int index = -1;
    float highestInt = 0f;
    for (int i=0; i!=Value.length; i++){
      if (Value[i][1]>highestInt){
        index = i;
        highestInt = Value[i][1];
      }
    }
    return index;
  }
  
  /**
   * estimates the noise of a chromatogram based on the dynamic noise level algorithm of Xu H. and Freitas M.A. BMC Bioinformatics 2010
   * @return the noise level
   */
  public float estimateChromatogramNoise(){
    float noise = 0f;
    Vector<Float> intensities = new Vector<Float>();
    for( int i = 0; i < ScanCount; i ++ ) {
    	if (this.Value[i][1]>0)
    		intensities.add(this.Value[i][1]);
    }
    Collections.sort(intensities);      
    float delta = 0f;
    float noiseImprovement = 1.7f;
    for (int i=1;i!=(intensities.size()/2);i++){
    	delta += intensities.get(i)-intensities.get(i-1);
    }
    delta = delta/(float)(intensities.size()/2);
    float[][] matrixValues = new float[intensities.size()-1][2];
    float[][] intensityValues = new float[intensities.size()-1][1];
    for (int j=0;j!=intensities.size()-1;j++){
    	matrixValues[j][0] = (j+1);
    	matrixValues[j][1] = 1f;
    	intensityValues[j][0] = intensities.get(j);
    }
      
      
    for (int i=(intensities.size()/2);i!=intensities.size();i++){
    	FloatMatrix matrix = new FloatMatrix(matrixValues);
    	matrix.m = i-1;
    	FloatMatrix transposed = matrix.transpose();
    	FloatMatrix product = transposed.times(matrix);
    	FloatMatrix intensityMatrix = new FloatMatrix(intensityValues);
    	intensityMatrix.m = i-1;
    	FloatMatrix result = product.inverse().times(transposed).times(intensityMatrix);
    	noise = (float)i*result.A[0][0]+result.A[1][0];
    	if (intensities.get(i)>(noise*noiseImprovement)){
    		break;
    	}
    }
    return noise;
  }
  
}
