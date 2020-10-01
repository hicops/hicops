# Welcome to the HiCOPS documentation

* TOC
{:_data/toc}

```liquid
<h2>{{ site.data.toc.docs_list_title }}</h2>
<ul>
   {% for item in site.data.toc.docs %}
      <li><a href="{{ item.url }}">{{ item.title }}</a></li>
   {% endfor %}
</ul>
```
