---
{title: sw_parstr, link: sw_parstr, summary: parses input string, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_parstr, folder: swfiles, mathjax: 'true'}

---
 
### Syntax
 
`parsed = sw_parstr(strIn)`
 
### Description
 
`parsed = sw_parstr(strIn)` parses input string containing a mathematical
expression, a linear combination of symbols with numerical prefactors.
The numerical symbols begin with `'S'`, `'M'` or `'P'` followed by two
letters from the set of `'xyz'` or `'P'` can be followed be a single
letter from the `'xyz'` set. For example a valid input string is
`'Sxx-0.33*Syy'`.
 
### Examples
 
Simple test:
```matlab
sw_parstr('Sxx + Syy')
```
*Output*
```
  struct with fields:
       type: {[2 1 1]  [2 2 2]}
    preFact: [1 1]
     string: 'Sxx + Syy'
```
 
 
### Input Arguments
 
`strIn`
: String that contains a linear combination of symbols, e.g
  `'0.1*Sxx+0.23*Syy'`.
 
### Output Arguments
 
`parsed`
: Struct type with the following fields:
  * `type`    Cell contains the same number of elements as in the sum. Each element
              is a vector as follows:
    * `type{idx}(1)`    Index of type of cross section:
      * 1 for `Sperp`,
      * 2 for `Sab`,
      * 3 for `Mab`,
      * 4 for `Pab`,
      * 5 for  `Pa`.
    * `type{idx}(2:3)`  Index of the component:
      * 1 for `x`,
      * 2 for `y`,
      * 3 for `z`.
  * `preFact` Vector contains the values of the prefactors in the sum.
  * `string`  Original input string.
 
### See Also
 
[sw_egrid](sw_egrid) \| [spinw.fitspec](spinw_fitspec)

{% include links.html %}
