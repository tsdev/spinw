---
{title: sw_notify( ), keywords: sample, summary: sends notification in OSX, sidebar: sw_sidebar,
  permalink: swfiles_sw_notify.html, folder: swfiles, mathjax: 'true'}

---
  sends notification in OSX
 
  SW_NOTIFY(message)
 
  SW_NOTIFY('option1',value1, ...)
 
 
  The function can be used to show a notification when a calculation is
  finished. Note that the notification only shows up if Matlab is not
  the active window.
 
  Install terminal-notifier:
  Using Ruby, you can install it through RubyGems:
  $ [sudo] gem install terminal-notifier
  You can also install it via Homebrew:
  $ brew install terminal-notifier
 
  Options:
 
