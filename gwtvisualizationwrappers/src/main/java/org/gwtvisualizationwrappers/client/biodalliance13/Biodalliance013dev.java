package org.gwtvisualizationwrappers.client.biodalliance13;

import java.util.List;

import com.google.gwt.core.client.JavaScriptObject;
import com.google.gwt.core.client.ScriptInjector;

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
	 * Construct and show a Biodalliance genome browser visualization.
	 * @param urlPrefix url prefix, where to find Biodalliance resources (ie //www.biodalliance.org/release-0.14/)
	 * @param containerId div element id to inject Biodalliance visualization
	 * @param initChr init chr
	 * @param initViewStart init view start
	 * @param initViewEnd init view end
	 * @param currentConfig genome browser configuration
	 * @param sources source tracks
	 */
	public void show(
			String urlPrefix,
			String containerId, 
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
		    ScriptInjector.fromString(BiodallianceClientBundle.INSTANCE.polyfillFetch().getText())
		        .setWindow(ScriptInjector.TOP_WINDOW)
		        .inject();
		    
		    _init013Dev();
		}

		JavaScriptObject config = createNewBiodallianceBrowserConfig(urlPrefix, containerId, initChr, initViewStart, initViewEnd, currentConfig);
		
		//add a source(s)
		for (BiodallianceSource source : sources) {
			if (BiodallianceSource.SourceType.BIGWIG.equals(source.getSourceType())) {
				addBigwigSource(config, source);	
			} else if (BiodallianceSource.SourceType.VCF.equals(source.getSourceType())) {
				addVcfSource(config, source);
			} else if (BiodallianceSource.SourceType.BED.equals(source.getSourceType())) {
				addBedSource(config, source);
			}
		}
		
		_initBrowser(config);
	}

	private static native void _init013Dev() /*-{
		$wnd.Browser013Dev = $wnd.Browser;
		if ($wnd.fetch) {
			//may have been polyfilled into $wnd by fetch.js (for browsers that do not support fetch)
			self.fetch = $wnd.fetch;
		}
		
		$wnd.resolverFunction = function(url) {
		   return fetch(url, {  
			  credentials: 'include'  //sending credentials with a fetch request (session cookie)
			}).then(function(resp) {
		       return resp.json();
		   }).then(function(rdata) {
		       return rdata.url;
		   });
		}
		
	}-*/;

	private JavaScriptObject createNewBiodallianceBrowserConfig(
			String prefix,
			String containerId,
			String initChr,
			int initViewStart,
			int initViewEnd,
			BiodallianceConfigInterface config
			) {
		
		return createNewBiodallianceBrowserConfig(prefix, containerId, initChr, initViewStart, initViewEnd, config.getTwoBitURI(),
				config.getBwgURI(), config.getStylesheetURI(), config.getTrixURI(), config.getTrixxURI(), config.getSpeciesName(), config.getTaxon(), config.getAuthName(),
				config.getVersion(), config.getUscsName());
	}
	
	private native JavaScriptObject createNewBiodallianceBrowserConfig(
			String urlPrefix,
			String containerId,
			String initChr,
			int initViewStart,
			int initViewEnd,
			String genomeFileURI,
			String gencodeBBFileURI,
			String gencodeXMLFileURI, //stylesheet
			String gencodeIndexFileURI,
			String gencodeIndexxFileURI,
			String coordSystemSpeciesName,
			int coordSystemTaxon,
			String coordSystemAuth,
			String coordSystemVersion,
			String coordSystemUcscName
			) /*-{
				
		var biodallianceBrowserConfig = {
				uiPrefix: urlPrefix,
				pageName: containerId,
				noPersist: true,
				chr: initChr, 
				viewStart:  initViewStart, 
				viewEnd:  initViewEnd, 
				cookieKey: coordSystemSpeciesName, 
				fullScreen: false,
				maxHeight: 15000,
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
					resolver: $wnd.resolverFunction}, 

					{name: 'GENCODE',
					bwgURI: gencodeBBFileURI,
					stylesheet_uri: gencodeXMLFileURI,
					collapseSuperGroups: true, 
					trixURI: gencodeIndexFileURI,
					trixxURI: gencodeIndexxFileURI,
					subtierMax:5,
					pinned:true,
					resolver: $wnd.resolverFunction
					}]
		};
		return biodallianceBrowserConfig;
	}-*/;
	
	private void addBigwigSource(
			JavaScriptObject biodallianceBrowserConfig,
			BiodallianceSource source) {
		addBigwigSource(biodallianceBrowserConfig, source.getSourceName(), source.getSourceURI(), source.getStyleType(), source.getStyleGlyphType(), source.getStyleColor(), source.getHeightPx());
	}
	private native void addBigwigSource(
			JavaScriptObject biodallianceBrowserConfig,
			String sourceName,
			String sourceBwgURI,
			String styleType, 
			String styleGlyphType,
			String styleColor,
			int heightPx
			) /*-{
		var newSource = {
	    	name: sourceName,
			collapseSuperGroups:true,
			bwgURI: sourceBwgURI,
			style: [{type : styleType,
					style: {glyph: styleGlyphType,
							COLOR1:styleColor,
							COLOR2:styleColor,
							COLOR3:styleColor,
							HEIGHT:heightPx}}],
			resolver: $wnd.resolverFunction
	    }
	    biodallianceBrowserConfig.sources.push(newSource);
	}-*/;
	
	private void addVcfSource(
			JavaScriptObject biodallianceBrowserConfig,
			BiodallianceSource source) {
		addVcfSource(biodallianceBrowserConfig, source.getSourceName(), source.getSourceURI(), source.getSourceIndexURI(), source.getStyleType(), source.getStyleGlyphType(), source.getStyleColor(), source.getHeightPx(), source.getSourceType().name().toLowerCase());
	}
	
	private native void addVcfSource(
			JavaScriptObject biodallianceBrowserConfig,
			String sourceName,
			String sourceURI,
			String sourceIndexURI,
			String styleType, 
			String styleGlyphType,
			String styleColor,
			int heightPx,
			String sourcePayload
			) /*-{
		
	    var newSource = {
			name: sourceName,
			collapseSuperGroups:true,
			uri: sourceURI,
			indexURI: sourceIndexURI,
			payload: sourcePayload,
			tier_type: 'tabix',
			style: [ 
				{type: 'default',
				   style: {
				      glyph: 'PLIMSOLL',
				      HEIGHT:heightPx,
				      STROKECOLOR: 'black',
				      FGCOLOR: styleColor,
				      BUMP: true}}],
			resolver: $wnd.resolverFunction
	    }
	    biodallianceBrowserConfig.sources.push(newSource);
	}-*/;
	
	private void addBedSource(
			JavaScriptObject biodallianceBrowserConfig,
			BiodallianceSource source) {
		addBedSource(biodallianceBrowserConfig, source.getSourceName(), source.getSourceURI(), source.getSourceIndexURI(), source.getStyleType(), source.getStyleGlyphType(), source.getStyleColor(), source.getHeightPx(), source.getSourceType().name().toLowerCase());
	}

	
	private native void addBedSource(
			JavaScriptObject biodallianceBrowserConfig,
			String sourceName,
			String sourceURI,
			String sourceIndexURI,
			String styleType, 
			String styleGlyphType,
			String styleColor,
			int heightPx,
			String sourcePayload
			) /*-{
		var newSource = {
			name: sourceName,
			collapseSuperGroups:true,
			uri: sourceURI,
			indexURI: sourceIndexURI,
			payload: sourcePayload,
			tier_type: 'tabix',
			style: [ 
				{type: 'default',
				   style: {
				      glyph: 'BOX',
				      HEIGHT:heightPx,
				      STROKECOLOR: 'black',
				      BGCOLOR: styleColor,
				      BUMP: true,
				      LABEL: true}}],
			resolver: $wnd.resolverFunction
	    }
	    biodallianceBrowserConfig.sources.push(newSource);
	}-*/;
	
}
