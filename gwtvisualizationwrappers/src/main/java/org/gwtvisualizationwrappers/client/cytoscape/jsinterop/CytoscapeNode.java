package org.gwtvisualizationwrappers.client.cytoscape.jsinterop;

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

import com.google.gwt.core.client.js.JsType;

@JsType
public class CytoscapeNode {
	/** Automatically treated as @JsProperty **/	
	// whether the element is selected (default false)
	public boolean selected;
	// whether the selection state is mutable (default true)
	public boolean selectable;
	// when locked a node's position is immutable (default false)
	public boolean locked;
	// whether the node can be grabbed and moved by the user
	public boolean grabbable;
	// a space separated list of css class names that the element has
	public String classes;
	// the model position of the node (optional on init, mandatory after)
	public CytoscapeNodeXY position;
	// element data
	public CytoscapeNodeData data;
	
}
