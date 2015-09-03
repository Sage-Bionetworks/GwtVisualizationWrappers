package org.gwtvisualizationwrappers.client.biodalliance;

import java.util.Date;
import java.util.List;

import com.google.gwt.core.client.GWT;
import com.google.gwt.core.client.JavaScriptObject;
import com.google.gwt.core.client.JsonUtils;
import com.google.gwt.core.client.Scheduler;
import com.google.gwt.core.client.ScriptInjector;
import com.google.gwt.dom.client.Element;

/*
 * #%L
 * GwtVisualizationWrapper
 * %%
 * Copyright (C) 2015 GwtCytoscapeJs
 * %%
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * #L%
 */

public class Biodalliance013dev {
    
	/**
     * Check to see if Biodalliance JS version has been loaded already.
     * 
     * @return true if Biodalliance is loaded, false otherwise.
     */
    private native boolean isBiodalliance013DevLoaded() /*-{
    	return (typeof $wnd['Browser013Dev'] !== 'undefined');
    }-*/;
    
    
	private static native void _initBrowser(JavaScriptObject biodallianceBrowserConfig) /*-{
		$wnd.Browser = $wnd.Browser013Dev;
		var newBrowser = new $wnd.Browser(biodallianceBrowserConfig);
	}-*/;
	
	/**
	 * Construct and show a cytoscape graph.
	 * 
	 * @param containerId Element ID to put the graph into.
	 * @param cytoscapeGraphJson Exported JSON from Cytoscape.
	 * http://wiki.cytoscape.org/Cytoscape_3/UserManual#Cytoscape_3.2BAC8-UserManual.2BAC8-CytoscapeJs.Data_Exchange_between_Cytoscape_and_Cytoscape.js

{
	elements:{
		nodes:[],
		edges:[]
	}
	style:[{
			selector: 'node',
			style: {
			}
	}]
}
	 */
	/**
	 * 
	 * @param containerId div element id to inject Biodalliance visualization
	 * @param initChr init chr
	 * @param initViewStart init view start
	 * @param initViewEnd init view end
	 * @param currentConfig genome browser configuration
	 * @param sources source tracks
	 */
	public void show(String containerId, 
			String initChr,
			int initViewStart,
			int initViewEnd,
			BiodallianceConfigInterface currentConfig,
			List<BiodallianceSource> sources) {
		//lazy load the cytoscape.js source
		if (!isBiodalliance013DevLoaded()) {
		    ScriptInjector.fromString(BiodallianceClientBundle.INSTANCE.biodalliance0_13dev().getText())
		        .setWindow(ScriptInjector.TOP_WINDOW)
		        .inject();
		    _init013Dev();
		}

		JavaScriptObject config = createNewBiodallianceBrowserConfig(containerId, initChr, initViewStart, initViewEnd, currentConfig);
		
		//add a source(s)
		for (BiodallianceSource source : sources) {
			if (BiodallianceSource.SourceType.BIGWIG.equals(source.getSourceType())) {
				addBigwigSource(config, source);	
			} else if (BiodallianceSource.SourceType.VCF.equals(source.getSourceType())) {
				addVCFSource(config, source);
			}
		}
		
		_initBrowser(config);
	}

	private static native void _init013Dev() /*-{
		$wnd.Browser013Dev = $wnd.Browser;
	}-*/;

	
	

	private JavaScriptObject createNewBiodallianceBrowserConfig(
			String containerId,
			String initChr,
			int initViewStart,
			int initViewEnd,
			BiodallianceConfigInterface config
			) {
		
		return createNewBiodallianceBrowserConfig(containerId, initChr, initViewStart, initViewEnd, config.getTwoBitURI(),
				config.getBwgURI(), config.getStylesheetURI(), config.getTrixURI(), config.getSpeciesName(), config.getTaxon(), config.getAuthName(),
				config.getVersion(), config.getUscsName());
	}
	
	private native JavaScriptObject createNewBiodallianceBrowserConfig(
			String containerId,
			String initChr,
			int initViewStart,
			int initViewEnd,
			String genomeFileURI,
			String gencodeBBFileURI,
			String gencodeXMLFileURI, //stylesheet
			String gencodeIndexFileURI,
			String coordSystemSpeciesName,
			int coordSystemTaxon,
			String coordSystemAuth,
			String coordSystemVersion,
			String coordSystemUcscName
			) /*-{
		var resolverFunction = function(url) {
		   return fetch(url, {  
			  credentials: 'include'  //sending credentials with a fetch request (session cookie)
			}).then(function(resp) {
		       return resp.json();
		   }).then(function(rdata) {
		       return rdata.url;
		   });
		}
				
		var biodallianceBrowserConfig = {
				pageName: containerId,
				noPersist: true,
				chr: initChr, 
				viewStart:  initViewStart, 
				viewEnd:  initViewEnd, 
				cookieKey: coordSystemSpeciesName, 
				fullScreen: false,
				coordSystem: { 
					speciesName: coordSystemSpeciesName, 
					taxon: coordSystemTaxon, 
					auth: coordSystemAuth, 
					version: coordSystemVersion, 
					ucscName: coordSystemUcscName},
				baseColors: {
	                 'A': 'black',
	                 'C': 'black',
	                 'G': 'black',
	                 'T': 'black',
	                 '-': 'black', //deletion
	                 'I': 'red'    //insertion
	            },
				sources: [{name: 'Genome',
					twoBitURI: genomeFileURI,
					tier_type: 'sequence',
					provides_entrypoints: true,
					pinned: true,
					resolver: resolverFunction}, 

					{name: 'GENCODE',
					bwgURI: gencodeBBFileURI,
					stylesheet_uri: gencodeXMLFileURI,
					collapseSuperGroups: true, 
					trixURI: gencodeIndexFileURI,
					subtierMax:5,
					pinned:true,
					resolver: resolverFunction
					}]
		};
		return biodallianceBrowserConfig;
	}-*/;
	
	private void addBigwigSource(
			JavaScriptObject biodallianceBrowserConfig,
			BiodallianceSource source) {
		addBigwigSource(biodallianceBrowserConfig, source.getSourceName(), source.getSourceURI(), source.getStyleType(), source.getStyleGlyphType(), source.getStyleColor(), source.getTrackHeightPx());
	}
	private native void addBigwigSource(
			JavaScriptObject biodallianceBrowserConfig,
			String sourceName,
			String sourceBwgURI,
			String styleType, 
			String styleGlyphType,
			String styleColor,
			int trackHeightPx
			) /*-{
		var resolverFunction = function(url) {
		   return fetch(url, {  
			  credentials: 'include'  //sending credentials with a fetch request (session cookie)
			}).then(function(resp) {
		       return resp.json();
		   }).then(function(rdata) {
		       return rdata.url;
		   });
		}
	    var newSource = {
	    	name: sourceName,
			collapseSuperGroups:true,
			bwgURI: sourceBwgURI,
			style: [{type : styleType,
					style: {glyph: styleGlyphType,
							COLOR1:styleColor,
							COLOR2:styleColor,
							COLOR3:styleColor,
							HEIGHT:trackHeightPx}}],
			collapseSuperGroups: true, 
			resolver: resolverFunction
	    }
	    biodallianceBrowserConfig.sources.push(newSource);
	}-*/;
	
	private void addVCFSource(
			JavaScriptObject biodallianceBrowserConfig,
			BiodallianceSource source) {
		addVCFSource(biodallianceBrowserConfig, source.getSourceName(), source.getSourceURI(), source.getSourceIndexURI(), source.getStyleType(), source.getStyleGlyphType(), source.getStyleColor(), source.getTrackHeightPx());
	}
	
	private native void addVCFSource(
			JavaScriptObject biodallianceBrowserConfig,
			String sourceName,
			String sourceURI,
			String sourceIndexURI,
			String styleType, 
			String styleGlyphType,
			String styleColor,
			int trackHeightPx
			) /*-{
		var resolverFunction = function(url) {
		   return fetch(url, {  
			  credentials: 'include'  //sending credentials with a fetch request (session cookie)
			}).then(function(resp) {
		       return resp.json();
		   }).then(function(rdata) {
		       return rdata.url;
		   });
		}
	    var newSource = {
	    	name: sourceName,
			collapseSuperGroups:true,
			uri: sourceURI,
			indexURI: sourceIndexURI,
			payload: 'vcf',
			tier_type: 'tabix',
			style: [{type : styleType,
					style: {glyph: styleGlyphType,
							COLOR1:styleColor,
							COLOR2:styleColor,
							COLOR3:styleColor,
							HEIGHT:trackHeightPx}}],
			collapseSuperGroups: true, 
			resolver: resolverFunction
	    }
	    biodallianceBrowserConfig.sources.push(newSource);
	}-*/;
	
}
