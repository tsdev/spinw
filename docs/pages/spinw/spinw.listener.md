---
{title: spinw.listener method, link: spinw.listener, summary: LISTENER  Add listener
    for event without binding the listener to the source object., keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_listener, folder: spinw, mathjax: 'true'}

---

### Syntax

`    object hsource.  if hsource is an array of source handles, the listener`

### Description

    is a function handle that is invoked when the event is triggered.
 
    el = LISTENER(hSource, PropName, 'Eventname', callback) adds a 
    listener for a property event.  Eventname must be one of the character 
    vectors 'PreGet', 'PostGet', 'PreSet', and 'PostSet'.  PropName must be 
    either a single property name or cell array of property names, or an
    array of one ore more meta.property objects.  The properties must 
    belong to the class of hSource.  If hSource is scalar, PropName can 
    include dynamic properties.
    
    For all forms, listener returns an event.listener.  To remove a
    listener, delete the object returned by listener.  For example,
    delete(el) calls the handle class delete method to remove the listener
    and delete it from the workspace.  Calling delete(el) on the listener
    object deletes the listener, which means the event no longer causes
    the callback function to execute. 
 
    LISTENER does not bind the listener's lifecycle to the object that is
    the source of the event.  Destroying the source object does not impact
    the lifecycle of the listener object.  A listener created with LISTENER
    must be destroyed independently of the source object.  Calling delete(el) 
    explicitly destroys the listener. Redefining or clearing the variable
    containing the listener can delete the listener if no other references
    to it exist.  To tie the lifecycle of the listener to the lifecycle
    of the source object, use addlistener.
 

### See Also

[also] \| [addlistener] \| [event.listener] \| [spinw](spinw) \| [notify] \| [delete] \| [meta.property] \| [events]
Help for spinw/listener is inherited from superclass HANDLE
    Reference page in Doc Center
       doc spinw/listener

{% include links.html %}
