/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.msn.hydroxy.parser;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Properties;

import at.tugraz.genome.lda.exception.HydroxylationEncodingException;

/**
 * this class stores the letter encoding for the number of hydroxy groups
 * @author Juergen Hartler
 *
 */
public class HydroxyEncoding extends Properties{

  private static final long serialVersionUID = 4234382085044433873L;
  /** lookup from the textual encoding to the number of hydroxylation sites*/
  private Hashtable<String,Short> encodingToNumber_;
  
  /**
   * constructor requires the path to the lookup file
   * @param path path to hydroxylationEncoding.txt
   * @throws IOException exception that is thrown when the file is not there
   * @throws HydroxylationEncodingException when there is something wrong with the file
   */
  public HydroxyEncoding(String path) throws IOException, HydroxylationEncodingException{
    FileInputStream in = new FileInputStream(path);
    this.load(in);
    in.close();
    encodingToNumber_ = new Hashtable<String,Short>();
    String keyString = "";
    try {
      for (Object key : this.keySet()) {
        keyString = (String)key;
        if (keyString==null || keyString.length()==0)
          continue;
        encodingToNumber_.put((String)get(key),Short.parseShort(keyString));
      }
    }catch(NumberFormatException nfx) {
      throw new HydroxylationEncodingException("The value \""+keyString+"\" is not permitted as hydroxylation encoding - only integer numbers are allowed as keys!");
    }
  }
  
  /**
   * returns the encoded String for the number of hydroxy groups
   * @param hydroxyNumber the number of hydroxy groups
   * @return the encoded String for the number of hydroxy groups
   * @throws HydroxylationEncodingException 
   */
  public String getEncodedPrefix(short hydroxyNumber) throws HydroxylationEncodingException{
    String stringNumber = String.valueOf(hydroxyNumber);
    if (getProperty(stringNumber)==null)
      throw new HydroxylationEncodingException("The demanded hydroxylation number \""+hydroxyNumber+"\" does not exist in the encoding file");
    return getProperty(stringNumber);
  }

  /**
   * returns the number of hydroxy groups for the encoded String
   * @param encoded the encoded String for the number of hydroxy groups
   * @return the number of hydroxy groups
   * @throws HydroxylationEncodingException 
   */
  public Short getHydroxyNumber(String encoded) throws HydroxylationEncodingException{
    if (!encodingToNumber_.containsKey(encoded))
      throw new HydroxylationEncodingException("The demanded hydroxylation key \""+encoded+"\" does not exist in the encoding file");
    return encodingToNumber_.get(encoded);
  }
  
}
