from bokeh.plotting import figure, push
from bokeh.embed import autoload_server
from bokeh.session import Session
from bokeh.document import Document

# alternative to these lines, bokeh.io.output_server(...)
document = Document()
session = Session()
session.use_doc('population_reveal')
#session.load_document(document)
#
#plot = figure()
#plot.circle([1,2], [3,4])
#
#document.add(push)
#push(session, document)
#
#script = autoload_server(plot, session)
#print script
