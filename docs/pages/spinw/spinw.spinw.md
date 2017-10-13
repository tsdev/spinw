---
{title: spinw.spinw method, link: spinw.spinw, summary: Operators and special characters.,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_spinw, folder: spinw, mathjax: 'true'}

---

### Syntax

`arithmetic operators.`

### Description

  uplus      - Unary plus                         +    
  minus      - Minus                              -    
  uminus     - Unary minus                        -    
  mtimes     - Matrix multiply                    *    
  times      - Array multiply                    .*    
  mpower     - Matrix power                       ^    
  power      - Array power                       .^    
  mldivide   - Backslash or left matrix divide    \    
  mrdivide   - Slash or right matrix divide       /    
  ldivide    - Left array divide                 .\    
  rdivide    - Right array divide                ./    
  idivide    - Integer division with rounding option.
  kron       - Kronecker tensor product   
 
Relational operators.
  eq         - Equal                             ==     
  ne         - Not equal                         ~=     
  lt         - Less than                          <      
  gt         - Greater than                       >      
  le         - Less than or equal                <=     
  ge         - Greater than or equal             >=     
 
Logical operators.
  relop      - Short-circuit logical AND         &&     
  relop      - Short-circuit logical OR          ||     
  and        - Element-wise logical AND           &      
  or         - Element-wise logical OR            |      
  not        - Logical NOT                        ~      
  punct      - Ignore function argument or output ~
  xor        - Logical EXCLUSIVE OR
  any        - True if any element of vector is nonzero
  all        - True if all elements of vector are nonzero
 
Special characters. 
  colon      - Colon                              : 
  paren      - Parentheses and subscripting      ( )              
  paren      - Brackets                          [ ]     
  paren      - Braces and subscripting           { }          
  punct      - Function handle creation           @
  punct      - Decimal point                      .      
  punct      - Structure field access             .      
  punct      - Parent directory                   ..     
  punct      - Continuation                       ...    
  punct      - Separator                          ,      
  punct      - Semicolon                          ;      
  punct      - Comment                            %      
  punct      - Invoke operating system command    !            
  punct      - Assignment                         =
  punct      - Quote                              '      
  punct      - Double quote                       "    
  transpose  - Transpose                         .'
  ctranspose - Complex conjugate transpose        ' 
  horzcat    - Horizontal concatenation          [,]     
  vertcat    - Vertical concatenation            [;]     
  subsasgn   - Subscripted assignment          ( ),{ },.   
  subsref    - Subscripted reference           ( ),{ },.   
  numArgumentsFromSubscript - Number of arguments for indexing methods
  subsindex  - Subscript index
  metaclass  - Metaclass for MATLAB class         ?
 
Bitwise operators.
  bitand     - Bit-wise AND.
  bitcmp     - Complement bits.
  bitor      - Bit-wise OR.
  bitxor     - Bit-wise XOR.
  bitset     - Set bit.
  bitget     - Get bit.
  bitshift   - Bit-wise shift.
 
Set operators.
  union      - Set union.
  unique     - Set unique.
  intersect  - Set intersection.
  setdiff    - Set difference.
  setxor     - Set exclusive-or.
  ismember   - True for set member.
 

### See Also

[arith] \| [relop] \| [slash] \| [function_handle]

{% include links.html %}
