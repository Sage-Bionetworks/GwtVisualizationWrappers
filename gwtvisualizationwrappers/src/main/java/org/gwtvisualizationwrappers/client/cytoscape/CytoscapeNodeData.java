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

import com.google.gwt.core.client.js.JsType;

@JsType
public class CytoscapeNodeData {
  // element data
  //mandatory for each element, assigned automatically on undefined
  public String id;
  //indicates the compound node parent id; not defined => no parent
  public String parent;
  public String source;
  public String target;
}
