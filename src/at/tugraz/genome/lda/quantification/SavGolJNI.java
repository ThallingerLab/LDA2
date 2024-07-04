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

/**
 * SavGolJNI.java
 *
 * Java Native Interface for communication with the GPU for calculating the
 * Savitzky-Golay-Filter of a chromatogram
 * 
 * Created on: Feb 08, 2017
 *     Author: Stefan Brandstaetter
 *    Version: Feb 17, 2017
 */

public class SavGolJNI {
  static {
    System.loadLibrary("SavGol");
  }

  // address of the struct in the library, that holds the addresses of the device values,
  // to handle multiple threads
  // http://stackoverflow.com/a/19209870/6914637
  public long address_;

  /**
   * Native function that checks if a CUDA capable device is present
   * @return true if there is a CUDA capable device
   */
  public native boolean cudaCapableDeviceNative();

  /**
   * Calls the native function that allocates the memory on the GPU and handles the address
   * @param mallocSize - the size of the allocated memory
   */
  public void initMalloc( int mallocSize )
  {
    this.address_ = this.initMallocNative( mallocSize );
  }

  /**
   * Native function that allocates memory on the graphics card
   * @param mallocSize - the size of the allocated memory
   * @return the address of the allocated memory on the GPU
   */
  public native long initMallocNative( int mallocSize );

  /**
   * Calls the native function that smoothes the raw data
   * @param values contains the retention time and the raw data
   * @param numberOfScans - the length of the raw data 
   * @param range - the time range that includes datapoints that should be included
   * @param repeats - how often the chromatogram should be smoothed
   * @param startSmoothScan - start position of the smoothing
   * @param stopSmoothScan - stop position of the smoothing
   * @return
   */
  public float[] Smooth( float[][] values, int numberOfScans,
                         float range, int repeats, int startSmoothScan,
                         int stopSmoothScan )
  {
    return this.SmoothNative( values, numberOfScans, range, repeats,
                              startSmoothScan, stopSmoothScan, this.address_ );
  }

  /**
   * Native function that smoothes the raw data
   * @param values contains the retention time and the raw data
   * @param numberOfScans - the length of the raw data 
   * @param range - the time range that includes datapoints that should be included
   * @param repeats - how often the chromatogram should be smoothed
   * @param startSmoothScan - start position of the smoothing
   * @param stopSmoothScan - stop position of the smoothing
   * @param address of the allocated memory on the GPU
   * @return
   */
  public native float[] SmoothNative( float[][] values, int numberOfScans,
                                 float range, int repeats, int startSmoothScan,
                                 int stopSmoothScan, long address );

  /**
   * Calls the native function that frees the allocated space on the graphics card
   */
  public void Frees() {
      this.FreesNative(this.address_);
  }
  
  /**
   * Native function that frees the allocated space on the graphics card
   * @param address of the allocated memory on the GPU
   */
  public native void FreesNative( long address );
}