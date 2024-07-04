/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.msn.utils;

import java.io.IOException;

import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;

/**
 * helper classes for fragmentation rule sets
 * @author Juergen Hartler
 *
 */
public class RulesUtils
{
  
  /**
   * returns the amount of chain types; int[0] acyl/alkyl/alkenyl; int[1] LCB; int[2] acyl chains; int[3] alkyl chains; int[4] alkenyl chains
   * @param dir the directory where the rules reside in
   * @param rule the name of the rule
   * @return amount of chain types; int[0] acyl/alkyl/alkenyl; int[1] LCB; int[2] acyl chains; int[3] alkyl chains; int[4] alkenyl chains
   * @throws SpectrummillParserException exception if there is something wrong about the elementconfig.xml, or an element is not there
   * @throws IOException exception if there is something wrong about the file
   * @throws NoRuleException thrown if the rules are not there
   * @throws RulesException specifies in detail which rule has been infringed
   */
  public static int[] getAmountOfChainsCategorized (String dir, String rule) throws RulesException, NoRuleException, IOException, SpectrummillParserException {
    int[] chainAmounts = new int[5];
    int totalChains = Integer.parseInt(RulesContainer.getAmountOfChains(rule,dir));
    chainAmounts[1] = Integer.parseInt(RulesContainer.getAmountOfLCBs(rule, dir));
    chainAmounts[0] = totalChains-chainAmounts[1]; 
    chainAmounts[3] = Integer.parseInt(RulesContainer.getAmountOfAlkylChains(rule, dir));
    chainAmounts[4] = Integer.parseInt(RulesContainer.getAmountOfAlkenylChains(rule, dir));
    chainAmounts[2] = totalChains-chainAmounts[1]-chainAmounts[3]-chainAmounts[4];    
    return chainAmounts;
  }
  
}
