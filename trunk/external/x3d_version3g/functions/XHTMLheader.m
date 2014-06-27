function [data,loc_body]=XHTMLheader(options)
    [data,loc_html]=XMLaddNode('html');
    data=XMLaddProperty('xmlns','http://www.w3.org/1999/xhtml',data);
    [data,loc_head]=XMLaddNode('head',data,loc_html+1);
    data=XMLaddNode('meta',data,loc_head+1);
    data=XMLaddProperty('http-equiv','Content-Type',data);
    data=XMLaddNode('title',data,loc_head+1);
    data=XMLaddString(options.title,data);
    data=XMLaddNode('link',data,loc_head+1);
    data=XMLaddProperty('rel','stylesheet',data);
    data=XMLaddProperty('type','text/css',data);
    %data=XMLaddProperty('href','x3dom/x3dom.css',data);
    data=XMLaddNode('script',data,loc_head+1);
    data=XMLaddProperty('type','text/javascript',data);
    data=XMLaddProperty('src','http://www.x3dom.org/x3dom/release/x3dom.js',data);
    data=XMLaddString(' ',data);
    [data,loc_body]=XMLaddNode('body',data,loc_html+1);
    data=XMLaddNode('h1',data,loc_body+1);
    data=XMLaddString(options.title,data);
    if(options.interactive)
        data=XMLaddNode('p',data,loc_body+1);
        data=XMLaddString('Description',data);
        [data,loc_form]=XMLaddNode('form',data,loc_body+1);
        data=XMLaddProperty('method','post',data);
        data=XMLaddProperty('action','',data);
        data=XMLaddNode('textarea',data,loc_form+1);
        data=XMLaddProperty('id','comments',data);
        data=XMLaddProperty('name','ncomments',data);
        data=XMLaddProperty('cols','60',data);
        data=XMLaddProperty('rows','3',data);
    end
    
    

