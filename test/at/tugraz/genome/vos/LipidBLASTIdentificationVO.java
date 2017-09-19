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

package at.tugraz.genome.vos;

import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidBLASTIdentificationVO
{
  
  private Vector<LipidBLASTDetectionVO> detections_;
  
  private boolean sameFAs_;
  private boolean samePositions_;
  private Vector<LipidBLASTDetectionVO> detectionsOutsideRt_;
  
  public LipidBLASTIdentificationVO(){
    detections_ = new Vector<LipidBLASTDetectionVO>();
    sameFAs_ = false;
    detectionsOutsideRt_ = new Vector<LipidBLASTDetectionVO>();
  }
  
  
  public void addDetection(LipidBLASTDetectionVO detectVO){
    detections_.add(detectVO);
    sameFAs_ = true;
    samePositions_ = true;
    String sn1_ = null;
    String sn2_ = null;
    String sn3_ = null;
    Hashtable<String,String> validFAs = new Hashtable<String,String>();
    int count = 0;
    for (LipidBLASTDetectionVO detection : detections_){
      if (!sameFAs_ && !samePositions_) break;
      if (count==0){
        validFAs.put(detection.getSn1(), detection.getSn1());
        sn1_ = detection.getSn1();
        validFAs.put(detection.getSn2(), detection.getSn2());
        sn2_ = detection.getSn2();
        validFAs.put(detection.getSn3(), detection.getSn3());
        sn3_ = detection.getSn3();
        count++;
      } else {
        if (!validFAs.containsKey(detection.getSn1())) sameFAs_ = false;
        if (!sn1_.equalsIgnoreCase(detection.getSn1())) samePositions_ = false;
        if (!validFAs.containsKey(detection.getSn2())) sameFAs_ = false;
        if (!sn2_.equalsIgnoreCase(detection.getSn2())) samePositions_ = false;
        if (!validFAs.containsKey(detection.getSn3())) sameFAs_ = false;
        if (!sn3_.equalsIgnoreCase(detection.getSn3())) samePositions_ = false;
      }
    }
  }
  
  public Vector<String> getIdentifications(int amountSns){
    return this.getIdentifications(this.detections_,amountSns);
  }

  public Vector<String> getOutsideRtDetections(int amountSns){
    return this.getIdentifications(this.detectionsOutsideRt_,amountSns);
  }

  private Vector<String> getIdentifications(Vector<LipidBLASTDetectionVO> detections, int amountSns){
    Hashtable<String,String> usedFACombis = new Hashtable<String,String>();
    for (LipidBLASTDetectionVO detection : detections){
      String combiName = createCombiName(detection, amountSns);
      usedFACombis.put(combiName, combiName);
    }
    return new Vector<String>(usedFACombis.keySet());
  }
  
  public double getHighestProbability(String ident, int amountSns){
    double probab = 0d;
    for (LipidBLASTDetectionVO detection : detections_){
      String combiName = createCombiName(detection, amountSns);
      if (combiName.equalsIgnoreCase(ident) && detection.getProbability()>probab){
        probab = detection.getProbability();
      }
    }
    return probab;
  }
  
  private String createCombiName(LipidBLASTDetectionVO detection, int amountSns){
    String combiName = detection.getSn1();
    if (amountSns>1) combiName+="/"+detection.getSn2();
    if (amountSns>2) combiName+="/"+detection.getSn3();
    return combiName;
  }
  
  public Vector<String> getRetentionTimes(String ident, int amountSns){
    LinkedHashMap<String,String> rts = new LinkedHashMap<String,String>();
    for (LipidBLASTDetectionVO detection : detections_){
      String combiName = createCombiName(detection, amountSns);
      if (combiName.equalsIgnoreCase(ident)) rts.put(detection.getRetentionTime(),detection.getRetentionTime());
    }
    return new Vector<String>(rts.keySet());
  }

  public boolean hasSameFAs()
  {
    return sameFAs_;
  }

  public boolean hasSamePositions()
  {
    return samePositions_;
  }
  
  public void filterAwayWrongRts(double rt, double tol){
    Vector<LipidBLASTDetectionVO> newDetections = new Vector<LipidBLASTDetectionVO>();
    detectionsOutsideRt_ = new Vector<LipidBLASTDetectionVO>();
    for (LipidBLASTDetectionVO detect : detections_){
      double hitRt = Double.parseDouble(detect.getRetentionTime());
      if ((rt-tol)<hitRt && hitRt<(rt+tol)){
        newDetections.add(detect);
      }else{
        detectionsOutsideRt_.add(detect);
      }
    }
    this.detections_ = newDetections;
  }
  
  public boolean hasHitsWithinRt(){
    return this.detections_.size()>0;
  }
  
  public boolean hasRemovedRts(){
    return this.detectionsOutsideRt_.size()>0;
  }
  
  public Vector<LipidBLASTDetectionVO> getDetections(){
    return this.detections_;
  }
}
