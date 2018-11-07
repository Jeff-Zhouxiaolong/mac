function objs = emd_findtag(str)
//  locate objects with specific tag
// Calling Sequence
// emd_findtag(STR)
// Description
// Locate objects with specific tag STR
//
// emd_findtag(OBJECT_HANDLES,STR)
// Restricts the search to objects listed in objhandles and their descendants.
//
// emd_findtag(...,'-depth',d)
// The depth argument d controls how many levels under the handles in objhandles
// are traversed. Specify d as inf to get the default behavior of all levels.
// Specify d as 0 to restrict to the objects listed in OBJECT_HANDLES.
//
// Rem: In order for this to work properly, the object's tag field must be a string
// containing keywords (or tags) separated by commas.
//
//
// See also
//  emd_addtag
//  emd_hastag
//  emd_rmtag
// Authors
//  H. Nahrstaedt - Aug 2010
// G.Rilling 12/2006
// gabriel.rilling@ens-lyon.fr



//if or(typeof(varargin(1))=="handle")
  //lobj = varargin(1);
  //tag = regexptranslate('escape',varargin(2));

  //objs = findobj(lobj,'-regexp','Tag',['(\W|^)',tag,'(\W|$)'],varargin{3:end});
  objs = findobj("user_data",str);
// else
//   tag = regexptranslate('escape',varargin{1});
//
//   objs = findobj('-regexp','Tag',['(\W|^)',tag,'(\W|$)'],varargin{2:end});
// end
endfunction
