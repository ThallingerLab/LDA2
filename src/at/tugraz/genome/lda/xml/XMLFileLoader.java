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

import org.w3c.dom.*;
import javax.xml.parsers.*;
import org.xml.sax.*;

import at.tugraz.genome.lda.LipidDataAnalyzer;

import java.io.*;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class XMLFileLoader {

  protected File fn_;
  protected Document doc= null;
  
  public XMLFileLoader(File fn){
    fn_ = fn;
  }

  public void parseXMLFile() throws SAXException {
    boolean validation       = true;
    boolean ignoreWhitespace = true;
    boolean ignoreComments   = true;
    boolean putCDATAIntoText = false;
    boolean createEntityRefs = false;
    if (fn_.exists()){
      // Step 1: create a DocumentBuilderFactory and configure it
      DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
      // Optional: set various configuration options
      dbf.setValidating(validation);
      dbf.setIgnoringComments(ignoreComments);
      dbf.setIgnoringElementContentWhitespace(ignoreWhitespace);
      dbf.setCoalescing(putCDATAIntoText);
      // The opposite of creating entity ref nodes is expanding them inline
      dbf.setExpandEntityReferences(!createEntityRefs);
//      Schema schema = dbf.getSchema();
      // change location to look for DTDs
      String fullPath = LipidDataAnalyzer.class.getResource("LipidDataAnalyzer.class").toString();
      String pathToDTD = fullPath.substring(0,fullPath.length()-"at/tugraz/genome/lda/LipidDataAnalyzer.class".length())+"dtds/";
      // Step 2: create a DocumentBuilder that satisfies the constraints
      // specified by the DocumentBuilderFactory
      DocumentBuilder db = null;
      try {
        db = dbf.newDocumentBuilder();
      }catch (ParserConfigurationException pce) {
        System.err.println(pce);
        javax.swing.JOptionPane.showConfirmDialog(null,
            "Error Initalising Parser! \n Can not start program",
            "Parser Error",
            javax.swing.JOptionPane.CLOSED_OPTION,
            javax.swing.JOptionPane.ERROR_MESSAGE);
        System.exit(1);
      }
      // Set an ErrorHandler before parsing
      db.setErrorHandler(new MyErrorHandler());
      doc = null;
      try {
        doc = db.parse(new BufferedInputStream(new FileInputStream(fn_)),pathToDTD);
        db = null;
/*      } catch (SAXException se) {
        System.err.println(se.getMessage());
*/
      } catch (IOException ioe) {
        System.err.println(ioe);
      }
    }else{
      doc = null;
    }
  }



  public Document getDocument(){
    return this.doc;
  }



  // Error handler to report errors and warnings
  private static class MyErrorHandler implements ErrorHandler{
    /** Error handler output goes here */
//    private PrintWriter out;

    public MyErrorHandler() {
    }

    /**
     * Returns a string describing parse exception details
     */
    private String getParseExceptionInfo(SAXParseException spe){
      String systemId = spe.getSystemId();
      if (systemId == null) {
          systemId = "null";
      }
      String info = "URI=" + systemId +
          " Line=" + spe.getLineNumber() +
          ": " + spe.getMessage();
      return info;
    }

    // The following methods are standard SAX ErrorHandler methods.
    // See SAX documentation for more info.

    public void warning(SAXParseException spe) throws SAXException {
    }

    public void error(SAXParseException spe) throws SAXException {
      String message = "Error: " + getParseExceptionInfo(spe);
      throw new SAXException(message);
    }

    public void fatalError(SAXParseException spe) throws SAXException {
      String message = "Fatal Error: " + getParseExceptionInfo(spe);
      throw new SAXException(message);
    }
  }
}
