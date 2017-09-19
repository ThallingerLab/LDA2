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

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * 
 * @author Juergen Hartler
 * this class stores information about overlapping peaks
 */
public class IsoOverlapVO
{
  
  private String analyte_;
  private Hashtable<Integer,Vector<Integer>> isoAffectedBy_;
  private Vector<IsoLocationSpaces> overlapingPeaks_;
  
  /**
   * 
   * @param name of the analyte to be checked
   * @param maxIsos amount of isotopes detected for the analyte 
   */
  public IsoOverlapVO(String name, int maxIsos){
    analyte_ = name;
    this.initHashes();
  }
  
  private void initHashes(){
    isoAffectedBy_ = new Hashtable<Integer,Vector<Integer>>();
    overlapingPeaks_ = new Vector<IsoLocationSpaces>();    
  }
  
  /**
   * 
   * @return is there an overlap
   */
  public boolean hasOverlap(){
    if (overlapingPeaks_.size()==0) return false;
    else return true;
  }
  
  /**
   * adds an overlapping peak to the VO
   * @param overlapLocation
   */
  public void addOverlapLocation(IsoLocationSpaces overlapLocation){
    if (!isTheSameParam(overlapLocation)){
      int locNr = overlapingPeaks_.size();
      Vector<Integer> isos = overlapLocation.getAffectedIsos();
      for (Integer iso : isos){
        Vector<Integer> overlapLocs = new Vector<Integer>();
        if (isoAffectedBy_.containsKey(iso)) overlapLocs = isoAffectedBy_.get(iso);
        overlapLocs.add(locNr);
        isoAffectedBy_.put(iso, overlapLocs);
      }
      overlapingPeaks_.add(overlapLocation);
    }
  }
  
  /**
   * this method checks if there is already the same overlapping analyte but with a different name
   * @param overlapLocation
   * @return
   */
  private boolean isTheSameParam(IsoLocationSpaces overlapLocation){
    boolean isTheSameThere = false;
    for (IsoLocationSpaces overlap : overlapingPeaks_){
      if (overlap.isTheSame(overlapLocation))  isTheSameThere = true;
    }
    return isTheSameThere;
  }
  
  public String getId(){
    return this.analyte_;
  }
  
//  public void printOverlappingAnalytes(){
//    for (Integer iso : isoAffectedBy_.keySet()){
//      Vector<Integer> lapNr = isoAffectedBy_.get(iso);
//      String overlap = iso+": ";
//      for (Integer lap : lapNr){
//        IsoLocationSpaces space = overlapingPeaks_.get(lap);
//        overlap += space.getId()+" - ";
//        Hashtable<Integer,Vector<String>> laps = space.getOverlapsOfIso(iso);
//        for (Integer probeNr : laps.keySet()){
//          overlap += probeNr+"_";
//          for (String id : laps.get(probeNr)){
//            overlap += id+"/";
//          }
//          overlap = overlap.substring(0,overlap.length()-1);
//        }
//      }
//      System.out.println("     "+overlap);
//    }
//  }
  
  public boolean isCorrectionPossible(Hashtable<String,LipidParameterSet> correctedParams,Hashtable<String,String> removed){
    for (IsoLocationSpaces overlap : overlapingPeaks_){
      if (!(correctedParams.containsKey(overlap.getId())||removed.containsKey(overlap.getId()))) return false;
    }
    return true;
  }
  
  public LipidParameterSet makeIsotopicCorrection(LipidParameterSet toCorrect){
    Vector<Vector<CgProbe>> probesToCorrect =  toCorrect.getIsotopicProbes();
    Vector<Vector<CgProbe>> correctedProbes = new Vector<Vector<CgProbe>>();
    float totalArea = 0f;
    for (int iso=0; iso!=probesToCorrect.size();iso++){
      Vector<CgProbe> corrected = new Vector<CgProbe>();
      if (isoAffectedBy_.containsKey(iso)){
        Vector<CgProbe> uncorrected = probesToCorrect.get(iso);
        Vector<Integer> requLocs =  isoAffectedBy_.get(iso);
        for (int probeNr=0; probeNr!=uncorrected.size(); probeNr++){
          CgProbe probe = correctProbe(uncorrected.get(probeNr),iso,probeNr,requLocs);
          if (probe.Area>0) corrected.add(probe);
        }
      }else{
        corrected = probesToCorrect.get(iso);
      }
      for (CgProbe probe : corrected) totalArea += probe.Area;
      if (corrected.size()>0)correctedProbes.add(corrected);
      else break;
    }
    toCorrect.setIsotopicProbes(correctedProbes);
    toCorrect.Area = totalArea;
    return toCorrect;
  }
  
  private CgProbe correctProbe(CgProbe toCorrect, int iso, int probeNr, Vector<Integer> requLocs){
    boolean hasToBeCorrected = false;
    for (Integer locNr : requLocs){
      IsoLocationSpaces overlap = overlapingPeaks_.get(locNr);
      if (overlap.getOverlapsOfIso(iso,probeNr).size()>0){
        hasToBeCorrected=true;
        break;
      }
    }
    if (hasToBeCorrected){
      float area = toCorrect.Area;
      for (Integer locNr : requLocs){
        IsoLocationSpaces overlap = overlapingPeaks_.get(locNr);
        area -= overlap.getOtherIsotopeArea(iso,probeNr);
      }
      toCorrect.Area = area;
      return toCorrect;
    }else return toCorrect;
  }
  
  public boolean reinitOverlapLocations(LipidParameterSet set, Hashtable<String,LipidParameterSet> correctedParams, int maxIso){
    Vector<String> overlapIds = new Vector<String>();
    Hashtable<Integer,Boolean> negative = new Hashtable<Integer,Boolean>();
    Hashtable<Integer,Vector<Double>> distris = new Hashtable<Integer,Vector<Double>>();
    int count = 0;
    for (IsoLocationSpaces loc : overlapingPeaks_){
      overlapIds.add(loc.getId());
      negative.put(count, loc.getNegative());
      distris.put(count, loc.getDistri());
      count++;
    }
    this.initHashes();
    boolean hasStillOverlap = false;
    for (int i=0; i!=overlapIds.size();i++){
      String id = overlapIds.get(i);
      if (correctedParams.containsKey(id)){
        LipidParameterSet overlapSet = correctedParams.get(id);
        IsoLocationSpaces over = new IsoLocationSpaces(id,maxIso);
        over.setParameterSet(overlapSet,distris.get(i),negative.get(i));
        over.determineOverlappingParts(set);
        if (over.hasOverlap()){
          addOverlapLocation(over);
          hasStillOverlap = true;
        }
      }
    }
    return hasStillOverlap;
  }
}
