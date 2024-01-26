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
public class AbsoluteQuantSettingsWholeReader extends XMLFileLoader
{

  private AbsoluteQuantSettingsPanel quantSettingsPanel_;
  
  public AbsoluteQuantSettingsWholeReader(String fileName,AbsoluteQuantSettingsPanel quantSettingsPanel) throws SAXException
  {
    super(new File(fileName));
    quantSettingsPanel_ = quantSettingsPanel;
  }

  
  public void parseXMLFile() throws SAXException {
    super.parseXMLFile();
    if (doc!=null){
      Element root = doc.getDocumentElement();
      if (root.getNodeName().equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_ROOT)){
        NodeList rootList = root.getChildNodes();
        for (int i=0;i!=rootList.getLength();i++){
          Node rootSubNode = rootList.item(i);
          if (rootSubNode.getNodeName().equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_EXPS_SET)){
            quantSettingsPanel_.cleanVolumeSettings();
            NodeList experimentNodes = rootSubNode.getChildNodes();
            for (int j=0; j!=experimentNodes.getLength(); j++){
              if (experimentNodes.item(j).getNodeName().equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_EXP_SET))
                parseExperimentSettings((Element)experimentNodes.item(j));
            }
          }
          if (rootSubNode.getNodeName().equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_CLASSES_SET)){
            quantSettingsPanel_.cleanClassSettings();
            NodeList classesNodes = rootSubNode.getChildNodes();
            for (int j=0; j!=classesNodes.getLength(); j++){
              if (classesNodes.item(j).getNodeName().equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_CLASS_SET))
                parseClassesSettings((Element)classesNodes.item(j));
            }
          }          
        }
      }
    }
  }
  
  private void parseExperimentSettings(Element experiment){
    String name = experiment.getAttribute(XMLConstants.WHOLE_ABS_SET_EXP_SET_ATTR_NAME);
    ExpVolumeSettingsPanel volSet = quantSettingsPanel_.getVolumeSettings().get(name);
    if (volSet!=null){
      NodeList nl = experiment.getChildNodes();
      for (int i=0; i!=nl.getLength();i++){
        if (nl.item(i).getNodeName().equalsIgnoreCase(XMLConstants.TAG_SETTING)){
          Element setting = (Element) nl.item(i);
          String settingType = setting.getAttribute(XMLConstants.ATTR_NAME);
          if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_PROBE_VOL)){
            setSetting(volSet.getProbeVolume(),setting);
          }else if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_END_VOL)){
            setSetting(volSet.getEndVolume(),setting);
          }else if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_SAMPLE_WEIGHT)){
            setSetting(volSet.getSampleWeight(),setting);     
          }else if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_PROTEIN_CONC)){
            setSetting(volSet.getProteinConc(),setting); 
          }else if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_LIPID_CONC)){
            setSetting(volSet.getNeutralConc(),setting); 
          }
        }
      }
    }
  }
  
  private void parseClassesSettings(Element classSet){
    String className = classSet.getAttribute(XMLConstants.WHOLE_ABS_SET_CLASS_NAME);
    LipidClassSettingsPanel classPanel = quantSettingsPanel_.getClassSettings().get(className);
    if (classPanel == null) return;
    if (classPanel.getChosenClass()!=className)
    {
    	classPanel = quantSettingsPanel_.getClassSettings().get(classPanel.getChosenClass());
    }
    if (classPanel != null && classPanel.areStandardsAvailable()){
      NodeList nl1 = classSet.getChildNodes();
      for (int i=0; i!=nl1.getLength(); i++){
        String nodeName = nl1.item(i).getNodeName();
        if (nodeName.equalsIgnoreCase(XMLConstants.TAG_SETTING)){
          Element setting = (Element) nl1.item(i);
          String settingType = setting.getAttribute(XMLConstants.ATTR_NAME);
          if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_DILUTION))
            classPanel.getGeneralSettingsPanel().getDilutionFactor().setText(setting.getAttribute(XMLConstants.ATTR_VALUE));
          else if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_ALL_EXPS_SAME)){
            String value = setting.getAttribute(XMLConstants.ATTR_VALUE);
            if (value!=null && value.equalsIgnoreCase("true"))
              classPanel.getAllSettingsSame().setSelected(true);
            else
              classPanel.getAllSettingsSame().setSelected(false);
            classPanel.changeVisibilitySettingsSame();
          }
          else if (settingType.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_ALL_STANDS_SAME)){
            String value = setting.getAttribute(XMLConstants.ATTR_VALUE);
            if (value!=null && value.equalsIgnoreCase("true"))
              classPanel.getAllStandardsSame().setSelected(true);
            else
              classPanel.getAllStandardsSame().setSelected(false);
            classPanel.changeVisibilityStandardsSame();
          }
        } else if (nodeName.equalsIgnoreCase(XMLConstants.TAG_ES_SETTING)){
          Element esSetting = (Element) nl1.item(i);
          StandardsSettingsPanel standSet = classPanel.getGeneralSettingsPanel();
          setStandardSettingsPanel(standSet, esSetting, false);
        } else if (nodeName.equalsIgnoreCase(XMLConstants.TAG_IS_SETTING)){
          Element isSetting = (Element) nl1.item(i);
          StandardsSettingsPanel standSet = classPanel.getGeneralSettingsPanel();
          setStandardSettingsPanel(standSet, isSetting, true);          
        } else if (nodeName.equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_EXPS_STAND_SET)){
          NodeList nl2 = nl1.item(i).getChildNodes();
          for (int j=0; j!=nl2.getLength(); j++){
            if (nl2.item(j).getNodeName().equalsIgnoreCase(XMLConstants.WHOLE_ABS_SET_EXP_STAND_SET)){
              Element expSetting = (Element)nl2.item(j);
              String expName = expSetting.getAttribute(XMLConstants.ATTR_NAME);
              StandardsSettingsPanel standSet = classPanel.getStandardsSettings().get(expName);
              if (standSet != null){
                NodeList nl3 = expSetting.getChildNodes();
                for (int k=0; k!= nl3.getLength(); k++){
                  String nl3NodeName = nl3.item(k).getNodeName();
                  if (nl3NodeName.equalsIgnoreCase(XMLConstants.TAG_ES_SETTING)){
                    Element esSetting = (Element) nl3.item(k);
                    setStandardSettingsPanel(standSet, esSetting, false);
                  } else if (nl3NodeName.equalsIgnoreCase(XMLConstants.TAG_IS_SETTING)){
                    Element isSetting = (Element) nl3.item(k);
                    setStandardSettingsPanel(standSet, isSetting, true);                              
                  } else if (nl3NodeName.equalsIgnoreCase(XMLConstants.TAG_SETTING)){
                    String value = ((Element) nl3.item(k)).getAttribute(XMLConstants.ATTR_VALUE);
                    standSet.getDilutionFactor().setText(value);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  private void setStandardSettingsPanel(StandardsSettingsPanel panel, Element setting, boolean isIS){
    PhysicalUnitInput volumeInp = null;
    PhysicalUnitInput concInp = null;
    String name = setting.getAttribute(XMLConstants.ATTR_NAME);
    if (isIS){
      if (name.equalsIgnoreCase(XMLConstants.IS_SETTINGS_VALUE_GENERAL)){
        volumeInp = panel.getGeneralIntVolume();
        concInp = panel.getGeneralIntConc();
      }else{
        volumeInp = panel.getIntVolumeSettings().get(name);
        concInp = panel.getIntConcSettings().get(name);
      }
    } else {
      if (name.equalsIgnoreCase(XMLConstants.ES_SETTINGS_VALUE_GENERAL)){
        volumeInp = panel.getGeneralExtVolume();
        concInp = panel.getGeneralExtConc();
      }else{
        volumeInp = panel.getExtVolumeSettings().get(name);
        concInp = panel.getExtConcSettings().get(name);
      }      
    }
    if (volumeInp!=null && concInp!=null){
      NodeList nl = setting.getChildNodes();
      for (int i=0; i!=nl.getLength(); i++){
        if (nl.item(i).getNodeName().equalsIgnoreCase(XMLConstants.TAG_SETTING)){
          Element set = (Element)nl.item(i);
          String type = set.getAttribute(XMLConstants.ATTR_NAME);
          if (type.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_INT_VOLUME)||
              type.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_EXT_VOLUME))
            setSetting(volumeInp,set);
          else if (type.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_INT_CONC)||
              type.equalsIgnoreCase(XMLConstants.WHOLE_ABS_NAME_EXT_CONC))
            setSetting(concInp,set);
        }
      }
    }
  }
  
  private void setSetting(PhysicalUnitInput inputMask, Element setting){
    String value = setting.getAttribute(XMLConstants.ATTR_VALUE);
    String magn = setting.getAttribute(XMLConstants.TAG_SETTING_ATTR_MAGNITUDE);
    inputMask.setInputValue(value, magn);
  }
}
