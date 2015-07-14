package org.gwtvisualizationwrappers.client.biodalliance;

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

public class Biodalliance013 {
    
	/**
     * Check to see if Biodalliance JS version has been loaded already.
     * 
     * @return true if Biodalliance is loaded, false otherwise.
     */
    private native boolean isBiodalliance013Loaded() /*-{
    	return (typeof $wnd['Browser013'] !== 'undefined');
    }-*/;
    
    
	private static native void _initBrowser(String containerId) /*-{
		var biodallianceBrowserConfig = {
			pageName: containerId,
			chr: '21', 
			viewStart:  33031597, 
			viewEnd:  33041570, 
			cookieKey: 'human', 
			fullScreen: true,
			coordSystem: { 
				speciesName: 'human', 
				taxon: 9606, 
				auth: 'NCBI', 
				version: '37', 
				ucscName: 'hg19'},
			baseColors: {
                 'A': 'black',
                 'C': 'black',
                 'G': 'black',
                 'T': 'black',
                 '-': 'black', //deletion
                 'I': 'red'    //insertion
            },
			sources: 
				[{name: 'Genome',
				twoBitURI: 'Portal/filehandle?entityId=syn4557603&preview=false&proxy=false&version=1',
				tier_type: 'sequence',
				provides_entrypoints: true,
				pinned: true}, 

				{name: 'GENCODE',
				bwgURI: 'Portal/filehandle?entityId=syn4557576&preview=false&proxy=false&version=1', //'human/gencode.bb',
				stylesheet_uri: 'Portal/filehandle?entityId=syn4557577&preview=false&proxy=false&version=1',//'human/gencode.xml',	
				collapseSuperGroups: true, 
				trixURI: 'Portal/filehandle?entityId=syn4557578&preview=false&proxy=false&version=1',//'human/geneIndex.ix',
				subtierMax:5,
				pinned:true}
				,
				{name: 'A2_i14.mkdup.coordsort.bw',
					collapseSuperGroups:true,
					bwgURI: 'Portal/filehandle?entityId=syn3928320&preview=false&proxy=false&version=1',//'case/A2_i14.mkdup.coordsort.bw',
					style: [{type : 'default',
							style: {glyph: 'HISTOGRAM',
									COLOR1:'red',
									COLOR2:'red',
									COLOR3:'red',
									HEIGHT:30}}]
				}]
			};
		
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
	public void show(String containerId) {
		//lazy load the cytoscape.js source
		if (!isBiodalliance013Loaded()) {
		    ScriptInjector.fromString(BiodallianceClientBundle.INSTANCE.biodalliance0_13().getText())
		        .setWindow(ScriptInjector.TOP_WINDOW)
		        .inject();
		    _init013();
		}

		_initBrowser(containerId);
	}

	private static native void _init013() /*-{
		$wnd.Browser013 = $wnd.Browser;
		$wnd.Browser = undefined;
	}-*/;


}
