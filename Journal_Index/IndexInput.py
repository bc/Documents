import pickle
import numpy as np
from termcolor import cprint,colored

class index_topic:

	def __init__(self,topic,volume,page):
		assert type(int(volume))==int, colored("Volume must be an integer.",'green',attrs=['bold'])
		assert type(int(page))==int, colored("Page must be an integer.",'green',attrs=['bold'])
		self.topic = topic
		self.volume = volume
		self.page = page
		self.index = {str(volume):[page]}

	def add_pages_to_volume(self, volume, page):
		if str(volume) not in self.index.keys():
			self.index[str(volume)] = [page]
		else:
			self.index[str(volume)].append(page)
	def print_topic(self):
		keys = np.sort([int(key) for key in self.index.keys()])
		output = "&\\text{" + self.topic +"}"
		for j in range(len(keys)):
			pages = np.sort(self.index[str(keys[j])])
			discontinuity_indicator = (np.diff(pages)!=1)*np.arange(1,len(pages),1)
			if len(pages)==1:
				pages_string = "p. " + str(pages[0])
			else:
				pages_string = "pp. " + str(pages[0])
				if discontinuity_indicator[0] != 0: 
					pages_string += ", " + str(pages[1])
					if discontinuity_indicator[-1] != 0:
						for i in range(2,len(discontinuity_indicator)):
							if discontinuity_indicator[i] != 0:
								pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
					else:
						for i in range(2,len(discontinuity_indicator-1)):
							if discontinuity_indicator[i] != 0:
								pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
						pages_string += "-" + str(pages[-1])
				else:
					if discontinuity_indicator[-1] != 0:
						for i in range(1,len(discontinuity_indicator)):
							if discontinuity_indicator[i] != 0:
								pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
					else:
						for i in range(1,len(discontinuity_indicator-1)):
							if discontinuity_indicator[i] != 0:
								pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
						pages_string += "-" + str(pages[-1])
			if j == 0:
				output +=  "&& vol. " + str(keys[j]) + ": " + pages_string + "\\\\" + "\n"
			else:
				output += "& && vol. " + str(keys[j]) + ": " + pages_string + "\\\\" + "\n"
		return(output) 
def print_reference(topic_string,index):
	keys = np.sort([int(key) for key in index[topic_string].index.keys()])
	if len(topic_string)%2==1: 
		left_adjustment,right_adjustment = int((len(topic_string)-1)/2),int((len(topic_string)+1)/2)
	else:
		left_adjustment,right_adjustment = int(len(topic_string)/2),int(len(topic_string)/2)
	output = "\n" + " "*10 + colored(' '*(20-left_adjustment) + topic_string + ' '*(20-right_adjustment),'white',attrs=['bold','underline']) + " "*10
	for j in range(len(keys)):
		pages = np.sort(index[topic_string].index[str(keys[j])])
		discontinuity_indicator = (np.diff(pages)!=1)*np.arange(1,len(pages),1)
		if len(pages)==1:
			pages_string = "p. " + str(pages[0])
		else:
			pages_string = "pp. " + str(pages[0])
			if discontinuity_indicator[0] != 0: 
				pages_string += ", " + str(pages[1])
				if discontinuity_indicator[-1] != 0:
					for i in range(2,len(discontinuity_indicator)):
						if discontinuity_indicator[i] != 0:
							pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
				else:
					for i in range(2,len(discontinuity_indicator-1)):
						if discontinuity_indicator[i] != 0:
							pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
					pages_string += "-" + str(pages[-1])
			else:
				if discontinuity_indicator[-1] != 0:
					for i in range(1,len(discontinuity_indicator)):
						if discontinuity_indicator[i] != 0:
							pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
				else:
					for i in range(1,len(discontinuity_indicator-1)):
						if discontinuity_indicator[i] != 0:
							pages_string += "-" + str(pages[i]) + ", " + str(pages[i+1])
					pages_string += "-" + str(pages[-1])
		if j == 0:
			if (5+len(str(keys[j]))+2+len(pages_string))%2==1:
				left_adjustment,right_adjustment = int(((5+len(str(keys[j]))+2+len(pages_string))-1)/2),int(((5+len(str(keys[j]))+2+len(pages_string))+1)/2)
			else:
				left_adjustment,right_adjustment = int((5+len(str(keys[j]))+2+len(pages_string))/2),int((5+len(str(keys[j]))+2+len(pages_string))/2)
			output +=  colored("\n" + " "*(30-left_adjustment) + "vol. " + str(keys[j]) \
								+ ": " + pages_string + " "*(30-right_adjustment) \
								+ "\n",
								'white')
		else:
			if (5+len(str(keys[j]))+2+len(pages_string))%2==1:
				left_adjustment,right_adjustment = int(((5+len(str(keys[j]))+2+len(pages_string))-1)/2),int(((5+len(str(keys[j]))+2+len(pages_string))+1)/2)
			else:
				left_adjustment,right_adjustment = int((5+len(str(keys[j]))+2+len(pages_string))/2),int((5+len(str(keys[j]))+2+len(pages_string))/2)

			output += colored(" "*(30-left_adjustment) + "vol. " + str(keys[j]) \
								+ ": " + pages_string + " "*(30-right_adjustment)\
								+ "\n", \
								'white')
	cprint(output)
def list_references(index):
	cprint('-'*60,'white',attrs=['bold'])
	cprint('Reference List\n','white',attrs=['bold','underline'])
	[cprint('- ' + key,'white') for key in sorted(index.keys())]
	print('')
def print_single_reference(index):
	another_reference = True
	while another_reference == True:
		cprint('-'*60,'white',attrs=['bold'])
		cprint('Print Reference','blue',attrs = ['bold','underline'])
		valid_response_1 = False
		while valid_response_1 == False:
			topic_string = input(colored('Topic: ','red',attrs = ['bold'])).capitalize()
			if topic_string in index.keys():
				print_reference(topic_string,index)
				valid_response_1 = True
			elif topic_string == 'Exit':
				break
			else:
				cprint('Topic ' + topic_string + ' not in Index.','blue',attrs=['bold'])
				valid_response_1 = False
		valid_response_3 = False
		while valid_response_3 == False:
			cprint('-'*60,'white',attrs=['bold'])
			next_action = input(colored("Print Another Reference? ([y],n): ",'red',attrs=['bold'])).capitalize()
			if next_action in ['Y', '']:
				another_reference = True
				valid_response_3 = True
			elif next_action == 'N':
				another_reference = False
				valid_response_3 = True
			else:
				cprint('Invalid Response.', 'blue', attrs =['bold'])
				valid_response_3 = False
def delete_topic(index):
	delete_another_topic = True
	while delete_another_topic == True:
		cprint('-'*60,'white',attrs=['bold'])
		cprint('Delete Reference','blue',attrs = ['bold','underline'])
		valid_response_1 = False
		while valid_response_1 == False:
			topic_to_be_deleted = input(colored('Topic: ','red',attrs = ['bold'])).capitalize()
			if topic_to_be_deleted in index.keys():
				valid_response_2 = False
				while valid_response_2 == False:
					safety_check = input(colored("Are you sure you want to delete '" \
										+ topic_to_be_deleted + "' in Index? ([y],n): ",\
										'red',attrs = ['bold'])).capitalize()
					if safety_check in ['Y', '']:
						del index[topic_to_be_deleted]
						valid_response_2 = True
					elif safety_check == 'N':
						valid_response_2 = True
					else:
						cprint('Invalid Response.', 'blue', attrs =['bold'])
						valid_response_2 = False
				valid_response_1 = True
			elif topic_to_be_deleted == 'Exit':
				break
			else:
				cprint('Topic ' + topic_to_be_deleted + ' not in Index.','blue',attrs=['bold'])
				valid_response_1 = False
		valid_response_3 = False
		while valid_response_3 == False:
			cprint('-'*60,'white',attrs=['bold'])
			next_action = input(colored("Delete Another Reference? ([y],n): ",'red',attrs=['bold'])).capitalize()
			if next_action in ['Y', '']:
				delete_another_topic = True
				valid_response_3 = True
			elif next_action == 'N':
				delete_another_topic = False
				valid_response_3 = True
			else:
				cprint('Invalid Response.', 'blue', attrs =['bold'])
				valid_response_3 = False
def delete_volume(index):
	remove_another_volume_reference = True
	while remove_another_volume_reference == True:
		cprint('-'*60,'white',attrs=['bold'])
		cprint('Remove Volume Reference','blue',attrs = ['bold','underline'])
		valid_response_1 = False
		while valid_response_1 == False:
			topic_string = input(colored('Topic: ','red',attrs = ['bold'])).capitalize()
			if topic_string in index.keys():
				if len(index[topic_string].index.keys())==1:
					cprint("Topic '" + topic_string + "' has only one volume reference. \nSwitch to Delete Topic.",\
							'blue', attrs = ['bold'])
					valid_response_1 = True
				else:
					valid_response_2 = False
					while valid_response_2 == False:
						volume_to_be_deleted = input(colored('Volume: ','red',attrs = ['bold']))
						if volume_to_be_deleted in index[topic_string].index.keys():
							valid_response_3 = False
							while valid_response_3 == False:
								safety_check = input(colored("Are you sure you want to delete vol. " + volume_to_be_deleted + \
													" from topic '" + topic_string + " in Index? ([y],n): ",\
													'red',attrs = ['bold'])).capitalize()
								if safety_check in ['Y', '']:
									del index[topic_string].index[volume_to_be_deleted]
									valid_response_3 = True
								elif safety_check == 'N':
									valid_response_3 = True
								else:
									cprint('Invalid Response.', 'blue', attrs =['bold'])
									valid_response_3 = False
							valid_response_2 = True
						else:
							cprint("Topic '" + topic_string + "' is not in vol. " + volume_to_be_deleted + '.','blue',attrs=['bold'])
							valid_response_2 = False
					valid_response_1 = True
			elif topic_string == 'Exit':
				break
			else:
				cprint('Topic ' + topic_string + ' not in Index.','blue',attrs=['bold'])
				valid_response_1 = False
		valid_response_3 = False
		while valid_response_3 == False:
			cprint('-'*60,'white',attrs=['bold'])
			next_action = input(colored("Remove Another Volume Reference? ([y],n): ",'red',attrs=['bold'])).capitalize()
			if next_action in ['Y', '']:
				remove_another_volume_reference = True
				valid_response_3 = True
			elif next_action == 'N':
				remove_another_volume_reference = False
				valid_response_3 = True
			else:
				cprint('Invalid Response.', 'blue', attrs =['bold'])
				valid_response_3 = False
def delete_page(index):
	remove_another_page_reference = True
	while remove_another_page_reference == True:
		cprint('-'*60,'white',attrs=['bold'])
		cprint('Remove Page Reference','blue',attrs = ['bold','underline'])
		valid_response_1 = False
		while valid_response_1 == False:
			topic_string = input(colored('Topic: ','red',attrs = ['bold'])).capitalize()
			if topic_string in index.keys():
				valid_response_2 = False
				while valid_response_2 == False:
					volume_string = input(colored('Volume: ','red', attrs = ['bold']))
					if volume_string in index[topic_string].index.keys():
						if len(index[topic_string].index[volume_string])==1:
							cprint("Topic '" + topic_string + "' in vol. " + volume_string + \
									" has only 1 reference. \nSwitch to either Delete Topic or " + \
									"Remove Volume Reference.",\
									'blue',attrs=['bold'])
						else:
							valid_response_3 = False
							while valid_response_3 == False:
								page_to_delete = int(input(colored('Page: ','red',attrs=['bold'])))
								if page_to_delete in index[topic_string].index[volume_string]:
									valid_response_4 = False
									while valid_response_4 == False:
										safety_check = input(colored("Are you sure you want to remove the index of '" + topic_string +\
																		"' in vol. " + volume_string + ", p. " \
																		+ str(page_to_delete) + "? ([y],n): ",\
																		'red', attrs = ['bold'])).capitalize()
										if safety_check in ['Y', '']:
											index_of_page = index[topic_string].index[volume_string].index(page_to_delete)
											del index[topic_string].index[volume_string][index_of_page]
											valid_response_4 = True
										elif safety_check == 'N':
											valid_response_4 = True
										else:
											cprint('Invalid Response.', 'blue', attrs =['bold'])
											valid_response_4 = False
									valid_response_3 = True
								else:
									cprint("Page " + str(page_to_delete) + " in vol. " + volume_string \
												+ " does not include topic '" + topic_string + "'.", 'blue', attrs = ['bold'])
									valid_response_3 = False
						valid_response_2 = True
					else:
						cprint("Topic '" + topic_string + "' is not in vol. " + volume_string + '.','blue',attrs=['bold'])
						valid_response_2 = False
				valid_response_1 = True
			else:
				cprint("Topic '" + topic_string + "' not in Index.",'blue',attrs=['bold'])
				valid_response_1 = False
		valid_response_5 = False
		while valid_response_5 == False:
			cprint('-'*60,'white',attrs=['bold'])
			next_action = input(colored("Remove Another Page Reference? ([y],n): ",'red',attrs=['bold'])).capitalize()
			if next_action in ['Y', '']:
				remove_another_page_reference = True
				valid_response_5 = True
			elif next_action == 'N':
				remove_another_page_reference = False
				valid_response_5 = True
			else:
				cprint('Invalid Response.', 'blue', attrs =['bold'])
				valid_response_5 = False
def options(index):
	exit_options = False
	while exit_options == False:
		cprint('-'*60,'white',attrs=['bold'])
		cprint('Options','blue',attrs=['bold','underline'])
		cprint('\n- List All References          (1)'+\
			   '\n- Print Reference              (2)'+\
			   '\n- Delete Reference             (3)'+\
			   '\n- Delete Volume Reference      (4)'+\
			   '\n- Delete Page Reference        (5)'+\
			   '\n- [Exit]\n','blue',\
			   attrs=['bold'])
		cprint('-'*60,'white',attrs=['bold'])
		action_string = colored('Desired Action: ','red',attrs = ['bold'])
		action = input(action_string).capitalize()
		if action in ['Exit', '']:
			exit_options = True
		elif action == '1':
			list_references(index)
			exit_options = True
		elif action == '2':
			print_single_reference(index)
			exit_options = True
		elif action == '3':
			delete_topic(index)
			exit_options = True
		elif action == '4':
			delete_volume(index)
			exit_options = True
		elif action == '5':
			delete_page(index)
			exit_options = True
		else:
			cprint('Invalid Response.','blue',attrs=['bold'])
			exit_options = False

index = pickle.load(open('journalindeces.pkl','rb'))

additional_entry = True
while additional_entry == True:
	cprint('-'*60,'white',attrs=['bold'])
	topic_string = colored('Topic: ', 'red', attrs = ['bold'])
	topic = input(topic_string).capitalize()
	if topic == "Exit": break
	if topic == "Options": 
		options(index)
	else:
		volume_string = colored('Volume: ', 'red', attrs = ['bold'])
		volume = input(volume_string).capitalize()
		if volume == "Exit": break
		volume = int(volume)
		new_page = False
		while new_page == False:
			page_string = colored('Page: ', 'red', attrs = ['bold'])
			page = input(page_string).capitalize()
			if page == "Exit": break
			page = int(page)
			if topic in index.keys():
				if str(volume) not in index[topic].index.keys():
					new_page = True
				elif page in index[topic].index[str(volume)]:
					cprint("Page " + str(page) + " in vol. " + str(volume) \
							+ " is already indexed for topic '" + topic + "'.", 'blue', attrs = ['bold'])
					new_page = False
				else:
					new_page = True
			else:
				new_page = True

		if topic not in index.keys():
			index[topic] = index_topic(topic,volume,page)
		else:
			index[topic].add_pages_to_volume(volume,page)

	valid_entry = False
	while valid_entry == False:
		cprint('-'*60,'white',attrs=['bold'])
		response_string = colored("Additional Entries? ([y]/n): ",'red',attrs=['bold'])
		response = input(response_string).capitalize()
		if response not in ['Y','','N','Options']:
			cprint('Invalid Response. ','blue',attrs = ['bold'])
			valid_entry=False
		elif response in ['Y','']:
			additional_entry=True
			valid_entry=True
		elif response == 'Options':
			options(index)
			valid_entry = False
		else:
			cprint('-'*60,'white',attrs=['bold'])
			additional_entry=False
			valid_entry=True

pickle.dump(index,open('journalindeces.pkl','wb'),pickle.HIGHEST_PROTOCOL)

latexfile = open("main.tex","w")
latexfile.write("\\documentclass[a4paper]{article} \n" + \
					"\\usepackage[english]{babel} \n" + \
					"\\usepackage[utf8x]{inputenc} \n" + \
					"\\usepackage[T1]{fontenc} \n" + \
					"\\usepackage{ragged2e} \n" + \
					"\\usepackage{amsmath} \n" + \
					"\\usepackage[a4paper,top=3cm,bottom=2cm,left=3cm,right=3cm,marginparwidth=1.75cm]{geometry} \n" + \
					"\\begin{document} \n" + \
					"\\section*{Index} \n" + \
					"\\begin{align*} \n")
for key in sorted(index.keys()):
	latexfile.write(index[key].print_topic())
latexfile.write("\\end{align*} \n" +\
				"\\end{document}")
latexfile.close()


