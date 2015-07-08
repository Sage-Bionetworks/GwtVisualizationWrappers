package org.gwtvisualizationwrappers.client.cytoscape;

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

import com.google.gwt.core.client.JavaScriptObject;
import com.google.gwt.core.client.JsArray;
import com.google.gwt.core.client.js.JsType;
import com.google.gwt.dom.client.Element;

@JsType
public class CytoscapeInitParams {
	public Element container;
	public CytoscapeElements elements;
	
//	style: [ /* ... */ ],
//	layout: { name: 'grid' /* , ... */ },
//	ready: function(evt){ /* ... */ },
	
	//initial viewport
	public int zoom;
	public CytoscapeNodeXY pan;
	
//	 minZoom: 1e-50,
//	 maxZoom: 1e50,
	public boolean zoomingEnabled;
	public boolean userZoomingEnabled;
	public boolean panningEnabled;
	public boolean userPanningEnabled;
	public boolean boxSelectionEnabled;
//	  selectionType: (isTouchDevice ? 'additive' : 'single'),
	public int touchTapThreshold;
	public int desktopTapThreshold;
	public boolean autolock;
	public boolean autoungrabify;
	public boolean autounselectify;
	// rendering options:
	public boolean headless;
	public boolean styleEnabled;
	public boolean hideEdgesOnViewport;
	public boolean hideLabelsOnViewport;
	public boolean textureOnViewport;
	public boolean motionBlur;
//	motionBlurOpacity: 0.2,
//	wheelSensitivity: 1,
//	pixelRatio: 1,
//	initrender: function(evt){ /* ... */ },
//	renderer: { /* ... */ }
}
