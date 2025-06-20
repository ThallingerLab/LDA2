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
import org.w3c.dom.Node;

import at.tugraz.genome.lda.swing.AbsoluteQuantSettingsPanel;
import at.tugraz.genome.lda.swing.ExpVolumeSettingsPanel;
import at.tugraz.genome.lda.swing.LipidClassSettingsPanel;
import at.tugraz.genome.lda.swing.PhysicalUnitInput;
import at.tugraz.genome.lda.swing.StandardsSettingsPanel;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AbsoluteQuantSettingsWholeWriter extends XMLFileWriter
{
  String fileName_;
  AbsoluteQuantSettingsPanel quantSettingsPanel_;
  
  public AbsoluteQuantSettingsWholeWriter (String fileName,AbsoluteQuantSettingsPanel quantSettingsPanel){
    super();
    this.fileName_ = fileName;
    this.quantSettingsPanel_ = quantSettingsPanel;
  }
  
  public void writeXMLFile() throws ParserConfigurationException,TransformerConfigurationException, TransformerException,IOException{
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    factory.setValidating(true);
    DOMImplementation implement = builder.getDOMImplementation();
    DocumentType docType = implement.createDocumentType(XMLConstants.WHOLE_DTD_NAME,"ID",XMLConstants.WHOLE_DTD_NAME+".dtd");
    Document newDoc = implement.createDocument("",XMLConstants.WHOLE_ABS_SET_ROOT,docType);
    Element root = newDoc.getDocumentElement();
    Element experimentSettings = newDoc.createElement(XMLConstants.WHOLE_ABS_SET_EXPS_SET);
    Hashtable<String,ExpVolumeSettingsPanel> volumeSettings = quantSettingsPanel_.getVolumeSettings();
    for (String expName : volumeSettings.keySet()){
      ExpVolumeSettingsPanel volSetPanel = volumeSettings.get(expName);
      Element expSetting = newDoc.createElement(XMLConstants.WHOLE_ABS_SET_EXP_SET);
      expSetting.setAttribute(XMLConstants.WHOLE_ABS_SET_EXP_SET_ATTR_NAME, expName);
      experimentSettings.appendChild(expSetting);
      this.setPhysicalValuesOfSettingsTag(newDoc, expSetting, XMLConstants.WHOLE_ABS_NAME_PROBE_VOL, volSetPanel.getProbeVolume());
      this.setPhysicalValuesOfSettingsTag(newDoc, expSetting, XMLConstants.WHOLE_ABS_NAME_END_VOL, volSetPanel.getEndVolume());
      this.setPhysicalValuesOfSettingsTag(newDoc, expSetting, XMLConstants.WHOLE_ABS_NAME_SAMPLE_WEIGHT, volSetPanel.getSampleWeight());
      this.setPhysicalValuesOfSettingsTag(newDoc, expSetting, XMLConstants.WHOLE_ABS_NAME_PROTEIN_CONC, volSetPanel.getProteinConc());
      this.setPhysicalValuesOfSettingsTag(newDoc, expSetting, XMLConstants.WHOLE_ABS_NAME_LIPID_CONC, volSetPanel.getNeutralConc());
    }
    root.appendChild(experimentSettings);
    Element classSettings = newDoc.createElement(XMLConstants.WHOLE_ABS_SET_CLASSES_SET);
    Hashtable<String,LipidClassSettingsPanel> classSets =  quantSettingsPanel_.getClassSettings();
    for (String className : classSets.keySet()){
      LipidClassSettingsPanel classPanel = classSets.get(className);
      if (classPanel.getAllSettingsSame() == null) continue; //TODO: in this case it is a class without available standards that may be defined by another class, consider writing this out.
      Element classSetting = newDoc.createElement(XMLConstants.WHOLE_ABS_SET_CLASS_SET);
      classSetting.setAttribute(XMLConstants.WHOLE_ABS_SET_CLASS_NAME,className);
      boolean allExpsSame = classPanel.getAllSettingsSame().isSelected();
      boolean allStandsSame = classPanel.getAllStandardsSame().isSelected();
      this.setPhysicalValuesOfSettingsTag(newDoc,classSetting, XMLConstants.WHOLE_ABS_NAME_ALL_EXPS_SAME,String.valueOf(allExpsSame),"");
      this.setPhysicalValuesOfSettingsTag(newDoc,classSetting, XMLConstants.WHOLE_ABS_NAME_ALL_STANDS_SAME,String.valueOf(allStandsSame),"");      
      if (allExpsSame){
        this.generateStandardsSettingsTag(newDoc, classSetting, classPanel.getGeneralSettingsPanel(),allStandsSame);
      }else{
        Hashtable<String,StandardsSettingsPanel> expSettings = classPanel.getStandardsSettings();
        Element expStandSets = newDoc.createElement(XMLConstants.WHOLE_ABS_SET_EXPS_STAND_SET);
        for (String expName : expSettings.keySet()){
          StandardsSettingsPanel expSetting = expSettings.get(expName);
          Element expStandSet = newDoc.createElement(XMLConstants.WHOLE_ABS_SET_EXP_STAND_SET);
          expStandSet.setAttribute(XMLConstants.ATTR_NAME, expName);
          this.generateStandardsSettingsTag(newDoc, expStandSet, expSetting, allStandsSame);
          expStandSets.appendChild(expStandSet);
        }
        classSetting.appendChild(expStandSets);
      }
      classSettings.appendChild(classSetting);
    }
    root.appendChild(classSettings);
    super.writeXMLFile(fileName_, newDoc);
  }
  
  private void generateStandardsSettingsTag(Document doc, Node parentNode, StandardsSettingsPanel settings,boolean allStandsSame){
    if (settings.getDilutionFactor().getText()!=null && settings.getDilutionFactor().getText().length()>0){
      this.setPhysicalValuesOfSettingsTag(doc,parentNode, XMLConstants.WHOLE_ABS_NAME_DILUTION,settings.getDilutionFactor().getText(),"");
    }
    if (allStandsSame){
      Element esSetting = doc.createElement(XMLConstants.TAG_ES_SETTING);
      esSetting.setAttribute(XMLConstants.ATTR_NAME, XMLConstants.ES_SETTINGS_VALUE_GENERAL);
      this.setPhysicalValuesOfSettingsTag(doc, esSetting, XMLConstants.WHOLE_ABS_NAME_EXT_VOLUME, settings.getGeneralExtVolume());
      this.setPhysicalValuesOfSettingsTag(doc, esSetting, XMLConstants.WHOLE_ABS_NAME_EXT_CONC, settings.getGeneralExtConc());      
      parentNode.appendChild(esSetting);
      
      Element isSetting = doc.createElement(XMLConstants.TAG_IS_SETTING);
      isSetting.setAttribute(XMLConstants.ATTR_NAME, XMLConstants.IS_SETTINGS_VALUE_GENERAL);
      this.setPhysicalValuesOfSettingsTag(doc, isSetting, XMLConstants.WHOLE_ABS_NAME_INT_VOLUME, settings.getGeneralIntVolume());
      this.setPhysicalValuesOfSettingsTag(doc, isSetting, XMLConstants.WHOLE_ABS_NAME_INT_CONC, settings.getGeneralIntConc());      
      parentNode.appendChild(isSetting);
    }else{
      Hashtable<String,PhysicalUnitInput> extVolSettings = settings.getExtVolumeSettings();
      for (String esName : extVolSettings.keySet()){
        Element esSetting = doc.createElement(XMLConstants.TAG_ES_SETTING);
        esSetting.setAttribute(XMLConstants.ATTR_NAME, esName);
        this.setPhysicalValuesOfSettingsTag(doc, esSetting, XMLConstants.WHOLE_ABS_NAME_EXT_VOLUME, extVolSettings.get(esName));
        this.setPhysicalValuesOfSettingsTag(doc, esSetting, XMLConstants.WHOLE_ABS_NAME_EXT_CONC, settings.getExtConcSettings().get(esName));
        parentNode.appendChild(esSetting);
      }
      Hashtable<String,PhysicalUnitInput> intVolSettings = settings.getIntVolumeSettings();
      for (String isName : intVolSettings.keySet()){
        Element esSetting = doc.createElement(XMLConstants.TAG_IS_SETTING);
        esSetting.setAttribute(XMLConstants.ATTR_NAME, isName);
        this.setPhysicalValuesOfSettingsTag(doc, esSetting, XMLConstants.WHOLE_ABS_NAME_INT_VOLUME, intVolSettings.get(isName));
        this.setPhysicalValuesOfSettingsTag(doc, esSetting, XMLConstants.WHOLE_ABS_NAME_INT_CONC, settings.getIntConcSettings().get(isName));
        parentNode.appendChild(esSetting);
      }
    }
  }
  
  private void setPhysicalValuesOfSettingsTag(Document doc, Node parentNode, String name, PhysicalUnitInput input){
    this.setPhysicalValuesOfSettingsTag(doc, parentNode, name, input.getValue().getText(), (String)input.getUnitMagnitude().getSelectedItem());
  }
  
  private void setPhysicalValuesOfSettingsTag(Document doc, Node parentNode, String name, String value, String magnitude){
    if (value!=null && value.length()>0){
      Element el = doc.createElement(XMLConstants.TAG_SETTING);
      el.setAttribute(XMLConstants.ATTR_NAME, name);
      el.setAttribute(XMLConstants.TAG_SETTING_ATTR_MAGNITUDE, magnitude);
      el.setAttribute(XMLConstants.ATTR_VALUE, value);
      parentNode.appendChild(el);
    }
  }
}
