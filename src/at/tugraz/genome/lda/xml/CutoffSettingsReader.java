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

package at.tugraz.genome.lda.xml;

import java.io.File;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.swing.CutoffSettingsPanel;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class CutoffSettingsReader extends XMLFileLoader
{
  
  CutoffSettingsPanel panel_;
  
  public CutoffSettingsReader(String fileName, CutoffSettingsPanel panel){
    super(new File(fileName));
    this.panel_ = panel;
  }
  
  public void parseXMLFile() throws SAXException {
    super.parseXMLFile();
    if (doc!=null){
      Element root = doc.getDocumentElement();
      if (root.getNodeName().equalsIgnoreCase(XMLConstants.CUTOFF_SET_ROOT)){
        NodeList rootList = root.getChildNodes();
        for (int i=0;i!=rootList.getLength();i++){
          Node rootSubNode = rootList.item(i);
          if (rootSubNode.getNodeName().equalsIgnoreCase(XMLConstants.CUTOFF_ISOTOPE)){
            String isotopeString = rootSubNode.getTextContent();
            if (isotopeString==null || isotopeString.length()==0) throw new SAXException("The content of "+XMLConstants.CUTOFF_ISOTOPE+" must not be empty!");
            try{
              int isotope = Integer.parseInt(isotopeString);
              if (isotope<0) throw new SAXException("The content of "+XMLConstants.CUTOFF_ISOTOPE+" must not be a negative integer!");
              panel_.setMaxIsotope(isotope);
            }catch(NumberFormatException nfx){
              throw new SAXException("The content of "+XMLConstants.CUTOFF_ISOTOPE+" must not be empty!");
            }catch(SettingsException sx){}
          }
          if (rootSubNode.getNodeName().equalsIgnoreCase(XMLConstants.CUTOFF_CUTOFF)){
            Element cutoff = (Element) rootSubNode;
            String lipidClass = cutoff.getAttribute(XMLConstants.ATTR_CLASS);
            String valueString = cutoff.getAttribute(XMLConstants.ATTR_VALUE);
            if (lipidClass==null || lipidClass.length()==0) throw new SAXException("The attribute "+XMLConstants.ATTR_CLASS+" of "+XMLConstants.CUTOFF_CUTOFF+" must not be empty!");
            if (valueString==null || valueString.length()==0) throw new SAXException("The attribute "+XMLConstants.ATTR_VALUE+" of "+XMLConstants.CUTOFF_CUTOFF+" must not be empty!");
            try{
              panel_.setCutoff(lipidClass,valueString);
            }catch(NumberFormatException nfx){
              throw new SAXException("The attribute "+XMLConstants.ATTR_VALUE+" of "+XMLConstants.CUTOFF_CUTOFF+" must be double format! The parsed value "+valueString+" is not!");
            }catch(AbsoluteSettingsInputException asx){
            }catch(SettingsException sx){
              throw new SAXException("The attribute "+XMLConstants.ATTR_VALUE+" of "+XMLConstants.CUTOFF_CUTOFF+" must not be negative! The parsed value "+valueString+" is negative!");
            }
            
          }
        }  
      }else{
        throw new SAXException("The document is not valid! The root element must be of type \""+XMLConstants.CUTOFF_SET_ROOT+"\"!");
      }
    }
  }
  
}
