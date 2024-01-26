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

import java.io.IOException;
import java.util.Hashtable;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;

import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.Element;

import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class CutoffSettingsWriter extends XMLFileWriter
{
  
  String fileName_;
  int maxIsotope_;
  Hashtable<String,String> cutoffs_;
  
  public CutoffSettingsWriter(String fileName, int maxIsotope, Hashtable<String,String> cutoffs){
    this.fileName_ = fileName;
    this.maxIsotope_ = maxIsotope;
    this.cutoffs_ = cutoffs;
  }
  
  public void writeXMLFile() throws ParserConfigurationException,TransformerConfigurationException, TransformerException,IOException, AbsoluteSettingsInputException{
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    factory.setValidating(true);
    DOMImplementation implement = builder.getDOMImplementation();
    DocumentType docType = implement.createDocumentType(XMLConstants.CUTOFF_DTD_NAME,"ID",XMLConstants.CUTOFF_DTD_NAME+".dtd");
    Document newDoc = implement.createDocument("",XMLConstants.CUTOFF_SET_ROOT,docType);
    Element root = newDoc.getDocumentElement();
    Element isotopeSetting = newDoc.createElement(XMLConstants.CUTOFF_ISOTOPE);
    isotopeSetting.appendChild(newDoc.createTextNode(String.valueOf(maxIsotope_)));
    root.appendChild(isotopeSetting);
    for (String className : cutoffs_.keySet()){
      Element cutoff = newDoc.createElement(XMLConstants.CUTOFF_CUTOFF);
      cutoff.setAttribute(XMLConstants.ATTR_CLASS, className);
      try{Double.parseDouble(cutoffs_.get(className));}catch(NumberFormatException ex){throw new AbsoluteSettingsInputException("The cut-off of the class "+className+" is not correct! "+ex.getMessage());}
      cutoff.setAttribute(XMLConstants.ATTR_VALUE, cutoffs_.get(className));
      root.appendChild(cutoff);
    }
    
    super.writeXMLFile(fileName_, newDoc);
  }  
}
