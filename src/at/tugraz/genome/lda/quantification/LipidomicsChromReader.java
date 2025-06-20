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

import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.quantification.CgProbe;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsChromReader extends ChromatogramReader
{
  private boolean useCuda_;
  protected SavGolJNI sav_gol_jni_;
  
  /**
   * This constructor is used on the cluster; here all of the information is definable
   * @param headerFilePath path to the header file
   * @param indexFilePath path to the index file
   * @param retentionTimeFilePath path to the retention time file
   * @param chromatogramFilePath path to the chromatogram file
   * @param sparseData - are there sparse time points in MS1 -> interpolation
   * @param chromSmoothRange if there is interpolation - the amount of points depend on the smooth range
   * @param useCuda if a CUDA capable device is installed
   * @throws CgException
   */
  public LipidomicsChromReader(String headerFilePath, String indexFilePath, String retentionTimeFilePath, String chromatogramFilePath,
      boolean sparseData, float chromSmoothRange, boolean useCuda) throws CgException{
    super(headerFilePath, indexFilePath, retentionTimeFilePath, chromatogramFilePath, sparseData, chromSmoothRange);
    this.useCuda_ = useCuda;
  }

  /** reads an m/z profile from the chrom file and smooths it */
  protected Vector<CgChromatogram> readProfiles(Vector<CgProbe> probes, float mzTolerance, float timeTolerance,float maxTimeDeviation,
      float mzSmoothRange, int smoothRepeats, int msLevel, SavGolJNI sav_gol_jni) throws CgException{
    sav_gol_jni_ = sav_gol_jni; 
    return super.readProfiles(probes, mzTolerance, timeTolerance, maxTimeDeviation, mzSmoothRange, smoothRepeats, msLevel);
  }

  /** smoothing the single profiles */
  protected Vector<CgChromatogram> smoothSingleProfiles(Vector<CgProbe> probes, Hashtable<Integer,CgChromatogram> singleProfiles, float mzSmoothRange,
      int smoothRepeats, Hashtable<Integer,Float> startTime, Hashtable<Integer,Float> stopTime){
    Vector<CgChromatogram> profiles = new Vector<CgChromatogram>(); 
    for (int i=0;i!=probes.size();i++){
  	  LipidomicsChromatogram cx = new LipidomicsChromatogram(singleProfiles.get(i));
      cx.isProfile_ = true;
      cx.Mz = probes.get(i).Peak;
      cx.LowerMzBand = startTime.get(i);
      cx.UpperMzBand = stopTime.get(i);
      if (useCuda_){
          cx.Smooth(mzSmoothRange, smoothRepeats, sav_gol_jni_);
        } else {
          cx.Smooth(mzSmoothRange, smoothRepeats);
        }
      profiles.add(cx);
    }    
    return profiles;
  }
 
  public Hashtable<Integer,Integer> getMsmsNrOfScans(){
    return msmsNrOfScans_;
  }
  
}
