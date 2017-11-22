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

package at.tugraz.genome.lda.analysis.vos;

import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class IsoLocationSpaces
{
  
  private String overlapAnalyte_;
  private int maxIsos_;
  
  private int nrOfPeaks_;
  // the first identifier is the peak number, the second the isotope number, the vector contains start and stop range
  private Hashtable<Integer,Hashtable<Integer,Vector<Float>>> timeRange_;
  private Hashtable<Integer,Hashtable<Integer,Vector<Float>>> mzRange_;
  private Hashtable<Integer,Float> zeroAreas_;
  
  private Vector<Vector<CgProbe>> probes_;
  
  private Hashtable<Integer,Hashtable<Integer,Vector<String>>> overlaps_;
  private Vector<Double> distri_;
  
  private boolean negative_;
  
  /**
   * The constructor stores the name and the maximum of isotopes of a potentially overlapping peak
   * 
   * @param name
   * @param maxIsotope
   */
  public IsoLocationSpaces(String name, int maxIsotope){
    overlapAnalyte_ = name;
    maxIsos_ = maxIsotope;
    overlaps_ = null;
  }
  
  /**
   * Here the ranges where an overlap can occur are set
   * 
   * @param set the lipid parameter set of the overlapping peak
   * @param negative does the overlapping peak has a negative distribution
   */
  public void setParameterSet(LipidParameterSet set, Vector<Double> distri, boolean negative){
    this.negative_ = negative;
    distri_ = distri;
    probes_ = set.getIsotopicProbes();
    nrOfPeaks_ = probes_.get(0).size();
    timeRange_ = new Hashtable<Integer,Hashtable<Integer,Vector<Float>>>();
    mzRange_ = new Hashtable<Integer,Hashtable<Integer,Vector<Float>>>();
    zeroAreas_ = new Hashtable<Integer,Float>();
    int count = 0;
    for (CgProbe probe : probes_.get(0)){
      Hashtable<Integer,Vector<Float>> timeRangeOfSinglePeak = new Hashtable<Integer,Vector<Float>>();
      Hashtable<Integer,Vector<Float>> mzRangeOfSinglePeak = new Hashtable<Integer,Vector<Float>>();
      float zeroTimeStart = probe.LowerValley;
      float zeroTimeEnd = probe.UpperValley;
      float zeroMzStart = probe.Mz-probe.LowerMzBand;
      float zeroMzEnd = probe.Mz+probe.UpperMzBand;
      if (probe instanceof Probe3D){
        Probe3D dProbe = (Probe3D)probe;
        zeroTimeStart = dProbe.getEllipseTimePosition()-dProbe.getEllipseTimeStretch();
        zeroTimeEnd = dProbe.getEllipseTimePosition()+dProbe.getEllipseTimeStretch();
        zeroMzStart = dProbe.getEllipseMzPosition()-dProbe.getEllipseMzStretch();
        zeroMzEnd = dProbe.getEllipseMzPosition()+dProbe.getEllipseMzStretch();
      }
      float negate = 1f;
      if (negative) negate = -1f;
      for (int i=0;i!=maxIsos_;i++){
        Vector<Float> timeRange = new Vector<Float>();
        Vector<Float> mzRange = new Vector<Float>();
        timeRange.add(zeroTimeStart);
        timeRange.add(zeroTimeEnd);
        mzRange.add(zeroMzStart+(((float)i)*LipidomicsConstants.getNeutronMass()*negate)/((float)set.getCharge()));
        mzRange.add(zeroMzEnd+(((float)i)*LipidomicsConstants.getNeutronMass()*negate)/((float)set.getCharge()));
        timeRangeOfSinglePeak.put(i, timeRange);
        mzRangeOfSinglePeak.put(i, mzRange);
      }
      timeRange_.put(count, timeRangeOfSinglePeak);
      mzRange_.put(count, mzRangeOfSinglePeak);
      float meanZeroArea = probe.Area;
      if (probes_.size()>1){
        float oneArea = -1f;
        if (probes_.get(0).size()==1 && probes_.get(1).size()==1){
          oneArea = probes_.get(1).get(0).Area;
        }else{
          //check which probe is in the range of the zero isotope - peak top must be covered by zero range 
          for (CgProbe probe1 : probes_.get(1)){
            if (probe.LowerMzBand<probe1.Peak && probe1.Peak < probe.UpperMzBand){
              oneArea = probe.Area;
              break;
            }
          }
        }
        if (oneArea>0){
          meanZeroArea = (meanZeroArea+(oneArea*((float)distri.get(0).doubleValue()))/((float)distri.get(1).doubleValue()))/2f;
        }
      }
      zeroAreas_.put(count, meanZeroArea);
      count++;
    }
  }
  
  /**
   * 
   * @param the peak to be checked for an overlap
   * @return true if the potentially overlapping peak overlaps
   */
  public boolean determineOverlappingParts(LipidParameterSet set){
    boolean overlap = false;
    overlaps_ = new Hashtable<Integer,Hashtable<Integer,Vector<String>>>();
    Vector<Vector<CgProbe>> probes = set.getIsotopicProbes();
    int iso=0;
    for (Vector<CgProbe> isoProbes : probes){
      int probeNr = 0;
      Hashtable<Integer,Vector<String>> overlapIso = new Hashtable<Integer,Vector<String>>();
      for (CgProbe probe : isoProbes){
        // the first is the isotope, the second is the probe number
        float timeStart = probe.LowerValley;
        float timeEnd = probe.UpperValley;
        float mzStart = probe.Mz-probe.LowerMzBand;
        float mzEnd = probe.Mz+probe.UpperMzBand;
        if (probe instanceof Probe3D){
          Probe3D dProbe = (Probe3D)probe;
          timeStart = dProbe.getEllipseTimePosition()-dProbe.getEllipseTimeStretch();
          timeEnd = dProbe.getEllipseTimePosition()+dProbe.getEllipseTimeStretch();
          mzStart = dProbe.getEllipseMzPosition()-dProbe.getEllipseMzStretch();
          mzEnd = dProbe.getEllipseMzPosition()+dProbe.getEllipseMzStretch();
        }
        Vector<String> overlaps = new Vector<String>();
        for (int i=0;i!=this.nrOfPeaks_;i++){
          Hashtable<Integer,Vector<Float>> timeRangeOfSinglePeak = timeRange_.get(i);
          Hashtable<Integer,Vector<Float>> mzRangeOfSinglePeak = mzRange_.get(i);
          Vector<String> isoOverlapsOfTheSamePeak = new Vector<String>();
          for (int j=0; j!=timeRangeOfSinglePeak.size();j++){
            Vector<Float> timeRange = timeRangeOfSinglePeak.get(j);
            Vector<Float> mzRange = mzRangeOfSinglePeak.get(j);
            if (isWithinRange(timeStart,timeEnd, timeRange.get(0),timeRange.get(1))&&isWithinRange(mzStart,mzEnd, mzRange.get(0),mzRange.get(1))){
              overlap = true;
              // the first is the probe, the second one the isotope
              String id_Overlap = i+"_"+j;
              isoOverlapsOfTheSamePeak.add(id_Overlap);
            }
          }
          if (LipidomicsConstants.useMostOverlappingIsotopeOnly() && isoOverlapsOfTheSamePeak.size()>0){
            String bestId = null;
            float highestOverlapPercentage = 0f;
            float range = mzEnd-mzStart;
            for (String id_Overlap : isoOverlapsOfTheSamePeak){
              int j = Integer.parseInt(id_Overlap.split("_")[1]);
              Vector<Float> mzRange = mzRangeOfSinglePeak.get(j);
              float start = mzRange.get(0);
              float stop = mzRange.get(1);
              if (mzStart>start) start = mzStart;
              if (stop>mzEnd) stop = mzEnd;
              float percentOverlap = (stop-start)/range;
              if (percentOverlap>highestOverlapPercentage){
                highestOverlapPercentage = percentOverlap;
                bestId = id_Overlap;
              }
              
            }
            isoOverlapsOfTheSamePeak = new Vector<String>();
            isoOverlapsOfTheSamePeak.add(bestId);
          }
          overlaps.addAll(isoOverlapsOfTheSamePeak);
        }
        if (overlaps.size()>0) overlapIso.put(probeNr, overlaps);
        probeNr++;
      }
      if (overlapIso.size()>0) overlaps_.put(iso, overlapIso);
      iso++;
    }
    return overlap;
  }
  
  /**
   * before this method, the method determineOverlappingParts has to be called
   * @return have these two peaks an overlap 
   */
  public boolean hasOverlap(){
    if (overlaps_==null) return false;
    for (int i=0;i!=maxIsos_;i++){
      if (this.overlaps_.containsKey(i)) return true;
    }
    return false;
  }
  
  /**
   * before this method, the method determineOverlappingParts has to be called
   * @return the isotopes that are affected by the overlap
   */
  public Vector<Integer> getAffectedIsos(){
    Vector<Integer> isos = new Vector<Integer>(overlaps_.keySet());
    Collections.sort(isos);
    return isos;
  }
  
  /**
   * determines if two probes are overlapping (this is a very simplistic approach)
   * @param start1
   * @param end1
   * @param start2
   * @param end2
   * @return
   */
  private boolean isWithinRange (float start1, float end1, float start2, float end2){
    boolean isInRange = false;
    if ((start1<start2&&end1>start2)||(start1>start2&&end2>start1)){
      isInRange = true;
    }
    return isInRange;
  }
  
  public Hashtable<Integer,Vector<String>> getOverlapsOfIso(int iso){
    return this.overlaps_.get(iso);
  }
  
  public Vector<String> getOverlapsOfIso(int iso, int probeNr){
    Vector<String> overlapInfo = new Vector<String>();
    Hashtable<Integer,Vector<String>> overlaps = this.getOverlapsOfIso(iso);
    if (overlaps.containsKey(probeNr)) overlapInfo = overlaps.get(probeNr);
    return overlapInfo;
  }
  
  /**
   * 
   * @return id of the overlapping analyte
   */
  public String getId(){
    return this.overlapAnalyte_;
  }
  
  /**
   * this method checks if there is already the same overlapping analyte but with a different name
   * only the zero isotope is checked, since this is the only relevant for the existence of an overlap
   * @param overlapLocation
   * @return
   */
  public boolean isTheSame(IsoLocationSpaces overlapLocation){
    //first, we check if we have the same amount of probes
    if (overlapLocation.probes_.size()==this.probes_.size()){
      if (overlapLocation.probes_.get(0).size()!=this.probes_.get(0).size()) return false;
    }else return false;
    // if the size is the same, we have to check the probes individually
    overlapLocation.probes_.get(0);
    Vector<CgProbe> probes1 = overlapLocation.probes_.get(0);
    Vector<CgProbe> probes2 = this.probes_.get(0);
    for (CgProbe probe1 : probes1){
      boolean found = false;
      for (CgProbe probe2 : probes2){
        if (checkTheSame(probe1,probe2)){
          found = true;
          break;
        }
      }
      if (!found) return false;
    }
    return true;
  }
  
  private boolean checkTheSame(CgProbe p1, CgProbe p2){
    if (!(p1.Area*0.99f<p2.Area&&p2.Area<p1.Area*1.01f)) return false;
    Probe3D p31 = null;
    Probe3D p32 = null;
    if (p1 instanceof Probe3D) p31 = (Probe3D)p1;
    if (p2 instanceof Probe3D) p32 = (Probe3D)p2;
    if (p31!=null&&p32!=null){
      if (!(p31.getEllipseMzPosition()*0.99f<p32.getEllipseMzPosition()&&p32.getEllipseMzPosition()<p31.getEllipseMzPosition()*1.01f)) return false;
      if (!(p31.getEllipseMzStretch()*0.99f<p32.getEllipseMzStretch()&&p32.getEllipseMzStretch()<p31.getEllipseMzStretch()*1.01f)) return false;
      if (!(p31.getEllipseTimePosition()*0.99f<p32.getEllipseTimePosition()&&p32.getEllipseTimePosition()<p31.getEllipseTimePosition()*1.01f)) return false;
      if (!(p31.getEllipseTimeStretch()*0.99f<p32.getEllipseTimeStretch()&&p32.getEllipseTimeStretch()<p31.getEllipseTimeStretch()*1.01f)) return false;
    } else if (p31==null && p32==null){
      if (!(p1.LowerValley*0.99f<p2.LowerValley&&p2.LowerValley<p1.LowerValley*1.01f)) return false;
      if (!(p1.UpperValley*0.99f<p2.UpperValley&&p2.UpperValley<p1.UpperValley*1.01f)) return false;
      if (!(p1.LowerMzBand*0.99f<p2.LowerMzBand&&p2.LowerMzBand<p1.LowerMzBand*1.01f)) return false;
      if (!(p1.UpperMzBand*0.99f<p2.UpperMzBand&&p2.UpperMzBand<p1.UpperMzBand*1.01f)) return false;
    } else return false;
    return true;
  }
  
  public boolean getNegative(){
    return this.negative_;
  }
  
  public float getOtherIsotopeArea(int iso, int probeNr){
    Vector<String> affectedPeaks = getOverlapsOfIso(iso,probeNr);
    float area=0f;
    for (String peakId : affectedPeaks){
      int probRef = Integer.parseInt(peakId.substring(0,peakId.indexOf("_")));
      int isoRef = Integer.parseInt(peakId.substring(peakId.indexOf("_")+1));
      //this if is for wrongly entered chemical formula!
      if (distri_.size()>isoRef)
        area += (zeroAreas_.get(probRef)*distri_.get(isoRef).floatValue())/distri_.get(0).floatValue();
    }
    return area;
  }
  
  public Vector<Double> getDistri(){
    return distri_;
  }
}
