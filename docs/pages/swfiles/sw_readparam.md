---
{title: sw_readparam, link: sw_readparam, summary: parse input arguments, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_readparam, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`parsed = sw_readparam(format,Name,Value)`
  
### Description
  
`parsed = sw_readparam(format,Name,Vale)` parses name-value pair
arguments. The parsing is controlled by the given `format` input. The
name-value pairs are converted into the parsed struct which has field
names identical to the given parameter names and corresponding values
taken from the input. `format` can also define required dimensionality of
a given value and default values for select parameters.
 
`sw_readparam` is used in most of the method functions of [spinw](spinw).
  
### Input Arguments
  
`format`
: A struct with the following fields:
  * `fname` Field names, $$n_{param}$$ strings in cell.
  * `size` Required dimensions of the corresponding value in a cell of
    $$n_{param}$$ vectors. Negative integer means dimension has to match with
    any other dimension which has the identical negative integer.
  * `defval` Cell of $$n_{param}$$ values, provides default values for
    missing parameters.
  * `soft` Cell of $$n_{param}$$ logical values, optional. If `soft(i)` is
    true, in case of missing parameter value $$i$$, no warning will be
    given.
 

{% include links.html %}
