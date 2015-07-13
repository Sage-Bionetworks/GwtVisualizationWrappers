package org.gwtvisualizationwrappers.client.cytoscape;

import com.google.gwt.core.client.GWT;
import com.google.gwt.core.client.JavaScriptObject;
import com.google.gwt.core.client.JsonUtils;
import com.google.gwt.core.client.Scheduler;
import com.google.gwt.core.client.ScriptInjector;
import com.google.gwt.dom.client.Element;

/*
 * #%L
 * GwtCytoscapeJs
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

public class CytoscapeGraph242 {
    
	/**
     * Check to see if Cytoscape JS version has been loaded already.
     * 
     * @return true if Cytoscape is loaded, false otherwise.
     */
    private native boolean isCytoscape242Loaded() /*-{
        return typeof $wnd['jQuery'].fn.cytoscape242 !== 'undefined'
    }-*/;
    
	private static native void _init242() /*-{
		$wnd.cytoscape242 = $wnd.cytoscape;
		$wnd.cytoscape = undefined;
	}-*/;
	
	private static native void _initGraph(String containerId, String cyjs, String styleJson) /*-{
		var containerElement = $doc.getElementById(containerId);
		function readyFunction() {
			console.log('Cytoscape graph ready');
		}
	
		var options = $wnd.createPlainObject(cyjs, styleJson, containerElement, readyFunction);
		
		$wnd.cytoscape = $wnd.cytoscape242;
		$wnd.cytoscape(options); 
		$wnd.cytoscape = undefined;
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
	public void show(String containerId, String cyjs, String styleJson) {
		//lazy load the cytoscape.js source
		if (!isCytoscape242Loaded()) {
		    ScriptInjector.fromString(CytoscapeJsClientBundle.INSTANCE.cytoscape2_4_2().getText())
		        .setWindow(ScriptInjector.TOP_WINDOW)
		        .inject();
		    _init242();
		}

		_initGraph(containerId, cyjs, styleJson);
	}



}
