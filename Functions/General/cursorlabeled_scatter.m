function main()
figure; hold on; axis([0 10 0 10]);

p1=scatter(2,2);
p2=scatter(4,5);
set(p1,'displayname','dogshit');
set(p2,'displayname','catpiss');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'Updatefcn',@datatip_name);

end

function output_txt = datatip_name(~,event_obj)
tar=get(event_obj,'Target');
traceid=get(tar,'DisplayName');
output_txt={['Gene: ',traceid]};
end